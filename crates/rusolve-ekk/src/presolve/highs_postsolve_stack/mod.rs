//! HiGHS source: presolve/HighsPostsolveStack.h (whole), presolve/HighsPostsolveStack.cpp (lines 1-82, 129-147)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! The stack of presolve reductions and their `undo` (postsolve). Split by
//! contiguous `.cpp` line range:
//!
//! | file | HighsPostsolveStack.cpp lines |
//! |---|---|
//! | `mod.rs` (this) | 20-82, 129-147 + every inline of the `.h` |
//! | `undo_a.rs` | 84-127, 149-576 (FreeColSubstitution .. DuplicateRow) |
//! | `undo_duplicate_col.rs` | 578-1012 (DuplicateColumn::undo, okMerge) |
//! | `undo_fix.rs` | 1014-1393 (DuplicateColumn::undoFix, transformToPresolvedSpace, SlackColSubstitution::undo) |
//!
//! # Rust mapping
//!
//! * The reduction structs are `pub` (private in C++ only to hide them from
//!   users). Their `undo` methods live in the impl files as
//!   `impl FreeColSubstitution { pub fn undo(...) }`.
//! * `HighsDataStack reductionValues` stores each struct with
//!   [`HighsDataStackItem`] (field by field, declaration order; see
//!   `util/highs_data_stack.rs`). `push(std::vector<Nonzero>)` is
//!   `push_vec`, `pop(std::vector<Nonzero>&)` is `pop_vec`.
//! * `const HighsMatrixSlice<F>&` arguments of the push templates are
//!   `impl IntoIterator<Item = HighsSliceNonzero>`; HPresolve passes
//!   `self.get_row_vector(row).iter(self)` etc., or `HighsEmptySlice`.
//! * `undoIterateBackwards` is an associated function (it uses no member).
//!   The `#ifdef DEBUG_EXTRA` branch is not compiled in the reference build:
//!   port the `#else` branch only.
//! * `undo(..., thread_safe)`: the C++ binds references
//!   `reductionValues_`, `colValues_`, `rowValues_` either to copies or to
//!   the members. Rust: `let (mut reduction_values_, mut col_values_, mut
//!   row_values_) = if thread_safe { (clone + reset_position, clone, clone) }
//!   else { self.reduction_values.reset_position(); (mem::take x3) };` and,
//!   when `!thread_safe`, move the three back into `self` before returning
//!   (the member state after the call is then identical to C++).
//! * `printf` diagnostics guarded by `report_col >= 0` become `print!`.
//! * `assert(1 == 0)` sentinels: `if cfg!(debug_assertions) { panic!(..) }`.

pub mod fme;
pub mod undo_a;
pub mod undo_duplicate_col;
pub mod undo_fix;
pub mod undo_fme;
pub mod undo_zero_cost;

pub use fme::*;

use crate::prelude::*;
use crate::util::highs_data_stack::impl_highs_data_stack_item;

/// HiGHS 1.16 (#3218 #3250): HighsPostsolveStack.h::HighsPostsolveStack::OrigRowType
/// (1.16.0 line 50). The `orig_row_type` vector is maintained in every mode; only the
/// `mip_presolve_row_append` paths read it.
#[repr(u8)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OrigRowType {
    #[default]
    kOriginal = 0,
    kCut = 1,
    kAppended = 2,
}

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::RowType
#[repr(u8)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RowType {
    #[default]
    kGeq = 0,
    kLeq = 1,
    kEq = 2,
}

impl HighsDataStackItem for RowType {
    const SIZE: usize = 1;
    fn write_bytes(&self, data: &mut Vec<u8>) {
        data.push(*self as u8);
    }
    fn read_bytes(data: &[u8], p: &mut usize) -> Self {
        let v = match data[*p] {
            0 => RowType::kGeq,
            1 => RowType::kLeq,
            _ => RowType::kEq,
        };
        *p += 1;
        v
    }
}

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::Nonzero
///
/// `Nonzero() = default` leaves the fields uninitialised; Rust zeroes them.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Nonzero {
    pub index: HighsInt,
    pub value: f64,
}

impl Nonzero {
    pub fn new(index_: HighsInt, value_: f64) -> Self {
        Self {
            index: index_,
            value: value_,
        }
    }
}
impl_highs_data_stack_item!(Nonzero {
    index: HighsInt,
    value: f64
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::LinearTransform
#[derive(Clone, Copy, Debug, Default)]
pub struct LinearTransform {
    pub scale: f64,
    pub constant: f64,
    pub col: HighsInt,
}
impl_highs_data_stack_item!(LinearTransform {
    scale: f64,
    constant: f64,
    col: HighsInt
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::FreeColSubstitution
#[derive(Clone, Copy, Debug, Default)]
pub struct FreeColSubstitution {
    pub rhs: f64,
    pub col_cost: f64,
    pub row: HighsInt,
    pub col: HighsInt,
    pub row_type: RowType,
}
impl_highs_data_stack_item!(FreeColSubstitution {
    rhs: f64,
    col_cost: f64,
    row: HighsInt,
    col: HighsInt,
    row_type: RowType
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::DoubletonEquation
#[derive(Clone, Copy, Debug, Default)]
pub struct DoubletonEquation {
    pub coef: f64,
    pub coef_subst: f64,
    pub rhs: f64,
    pub subst_lower: f64,
    pub subst_upper: f64,
    pub subst_cost: f64,
    pub row: HighsInt,
    pub col_subst: HighsInt,
    pub col: HighsInt,
    pub lower_tightened: bool,
    pub upper_tightened: bool,
    pub row_type: RowType,
}
impl_highs_data_stack_item!(DoubletonEquation {
    coef: f64,
    coef_subst: f64,
    rhs: f64,
    subst_lower: f64,
    subst_upper: f64,
    subst_cost: f64,
    row: HighsInt,
    col_subst: HighsInt,
    col: HighsInt,
    lower_tightened: bool,
    upper_tightened: bool,
    row_type: RowType
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::EqualityRowAddition
#[derive(Clone, Copy, Debug, Default)]
pub struct EqualityRowAddition {
    pub row: HighsInt,
    pub added_eq_row: HighsInt,
    pub eq_row_scale: f64,
}
impl_highs_data_stack_item!(EqualityRowAddition {
    row: HighsInt,
    added_eq_row: HighsInt,
    eq_row_scale: f64
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::EqualityRowAdditions
#[derive(Clone, Copy, Debug, Default)]
pub struct EqualityRowAdditions {
    pub added_eq_row: HighsInt,
}
impl_highs_data_stack_item!(EqualityRowAdditions {
    added_eq_row: HighsInt
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::SingletonRow
#[derive(Clone, Copy, Debug, Default)]
pub struct SingletonRow {
    pub coef: f64,
    pub row: HighsInt,
    pub col: HighsInt,
    pub col_lower_tightened: bool,
    pub col_upper_tightened: bool,
}
impl_highs_data_stack_item!(SingletonRow {
    coef: f64,
    row: HighsInt,
    col: HighsInt,
    col_lower_tightened: bool,
    col_upper_tightened: bool
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::FixedCol
#[derive(Clone, Copy, Debug, Default)]
pub struct FixedCol {
    pub fix_value: f64,
    pub col_cost: f64,
    pub col: HighsInt,
    pub fix_type: HighsBasisStatus,
}
impl_highs_data_stack_item!(FixedCol {
    fix_value: f64,
    col_cost: f64,
    col: HighsInt,
    fix_type: HighsBasisStatus
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::RedundantRow
#[derive(Clone, Copy, Debug, Default)]
pub struct RedundantRow {
    pub row: HighsInt,
}
impl_highs_data_stack_item!(RedundantRow { row: HighsInt });

/// HiGHS 1.16 (#3218): HighsPostsolveStack.h::HighsPostsolveStack::ImpliedEquation (lines 217-223)
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpliedEquation {
    pub row: HighsInt,
    pub at_lower: bool,
}
impl_highs_data_stack_item!(ImpliedEquation {
    row: HighsInt,
    at_lower: bool
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::ForcingRow
#[derive(Clone, Copy, Debug, Default)]
pub struct ForcingRow {
    pub side: f64,
    pub row: HighsInt,
    pub row_type: RowType,
}
impl_highs_data_stack_item!(ForcingRow {
    side: f64,
    row: HighsInt,
    row_type: RowType
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::ForcingColumn
#[derive(Clone, Copy, Debug, Default)]
pub struct ForcingColumn {
    pub col_cost: f64,
    pub col_bound: f64,
    pub col: HighsInt,
    pub at_infinite_upper: bool,
    pub col_integral: bool,
}
impl_highs_data_stack_item!(ForcingColumn {
    col_cost: f64,
    col_bound: f64,
    col: HighsInt,
    at_infinite_upper: bool,
    col_integral: bool
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::ForcingColumnRemovedRow
#[derive(Clone, Copy, Debug, Default)]
pub struct ForcingColumnRemovedRow {
    pub rhs: f64,
    pub row: HighsInt,
}
impl_highs_data_stack_item!(ForcingColumnRemovedRow {
    rhs: f64,
    row: HighsInt
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::DuplicateRow
#[derive(Clone, Copy, Debug, Default)]
pub struct DuplicateRow {
    pub duplicate_row_scale: f64,
    pub duplicate_row: HighsInt,
    pub row: HighsInt,
    pub row_lower_tightened: bool,
    pub row_upper_tightened: bool,
}
impl_highs_data_stack_item!(DuplicateRow {
    duplicate_row_scale: f64,
    duplicate_row: HighsInt,
    row: HighsInt,
    row_lower_tightened: bool,
    row_upper_tightened: bool
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::DuplicateColumn
#[derive(Clone, Copy, Debug, Default)]
pub struct DuplicateColumn {
    pub col_scale: f64,
    pub col_lower: f64,
    pub col_upper: f64,
    pub duplicate_col_lower: f64,
    pub duplicate_col_upper: f64,
    pub col: HighsInt,
    pub duplicate_col: HighsInt,
    pub col_integral: bool,
    pub duplicate_col_integral: bool,
}
impl_highs_data_stack_item!(DuplicateColumn {
    col_scale: f64,
    col_lower: f64,
    col_upper: f64,
    duplicate_col_lower: f64,
    duplicate_col_upper: f64,
    col: HighsInt,
    duplicate_col: HighsInt,
    col_integral: bool,
    duplicate_col_integral: bool
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::SlackColSubstitution
#[derive(Clone, Copy, Debug, Default)]
pub struct SlackColSubstitution {
    pub rhs: f64,
    pub row: HighsInt,
    pub col: HighsInt,
}
impl_highs_data_stack_item!(SlackColSubstitution {
    rhs: f64,
    row: HighsInt,
    col: HighsInt
});

/// HiGHS 1.16 (#2962): HighsPostsolveStack.h::ZeroObjSingletonContinuousCol (lines 286-299)
#[derive(Clone, Copy, Debug, Default)]
pub struct ZeroObjSingletonContinuousCol {
    pub orig_row_lower: f64,
    pub orig_row_upper: f64,
    pub new_row_lb: f64,
    pub new_row_ub: f64,
    pub lb: f64,
    pub ub: f64,
    pub coef: f64,
    pub col: HighsInt,
    pub row: HighsInt,
}
impl_highs_data_stack_item!(ZeroObjSingletonContinuousCol {
    orig_row_lower: f64,
    orig_row_upper: f64,
    new_row_lb: f64,
    new_row_ub: f64,
    lb: f64,
    ub: f64,
    coef: f64,
    col: HighsInt,
    row: HighsInt
});

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::ReductionType
#[repr(u8)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ReductionType {
    #[default]
    kLinearTransform = 0,
    kFreeColSubstitution,
    kDoubletonEquation,
    kEqualityRowAddition,
    kEqualityRowAdditions,
    kSingletonRow,
    kFixedCol,
    kRedundantRow,
    kForcingRow,
    kForcingColumn,
    kForcingColumnRemovedRow,
    kDuplicateRow,
    kDuplicateColumn,
    kSlackColSubstitution,
    // HiGHS 1.16 #2962: appended last so the existing discriminants are unchanged
    kZeroObjSingletonContinuousCol,
    // HiGHS 1.16 #3218: inserted after kRedundantRow in HighsPostsolveStack.h (line 312, 1.16.0);
    // appended last here for the same reason
    kImpliedEquation,
    // HiGHS #3255 (mip_fm_presolve): HighsPostsolveStack.h lines 318-319 (1.16.0)
    kFourierMotzkinBlock,
    kFourierMotzkinObjCol,
}

/// HiGHS: HighsPostsolveStack.cpp::computeRowStatus (lines 74-82, 9 lines; file-static)
pub(crate) fn compute_row_status(dual: f64, row_type: RowType) -> HighsBasisStatus {
    if row_type == RowType::kEq {
        if dual < 0.0 {
            HighsBasisStatus::kUpper
        } else {
            HighsBasisStatus::kLower
        }
    } else if row_type == RowType::kGeq {
        HighsBasisStatus::kLower
    } else {
        HighsBasisStatus::kUpper
    }
}

/// HiGHS: HighsPostsolveStack.cpp::computeStatus(double, HighsBasisStatus&, double) (lines 129-137, 9 lines; file-static)
pub(crate) fn compute_status_inout(
    dual: f64,
    status: &mut HighsBasisStatus,
    dual_feasibility_tolerance: f64,
) -> HighsBasisStatus {
    if dual > dual_feasibility_tolerance {
        *status = HighsBasisStatus::kLower;
    } else if dual < -dual_feasibility_tolerance {
        *status = HighsBasisStatus::kUpper;
    }

    *status
}

/// HiGHS: HighsPostsolveStack.cpp::computeStatus(double, double) (lines 139-147, 9 lines; file-static)
pub(crate) fn compute_status(dual: f64, dual_feasibility_tolerance: f64) -> HighsBasisStatus {
    if dual > dual_feasibility_tolerance {
        HighsBasisStatus::kLower
    } else if dual < -dual_feasibility_tolerance {
        HighsBasisStatus::kUpper
    } else {
        HighsBasisStatus::kBasic
    }
}

/// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack
///
/// C++ initialisers: `debug_prev_* = 0`, `origNumCol = -1`,
/// `origNumRow = -1`; everything else empty.
#[derive(Clone, Debug)]
pub struct HighsPostsolveStack {
    pub debug_prev_numreductions: usize,
    pub debug_prev_col_lower: f64,
    pub debug_prev_col_upper: f64,
    pub debug_prev_row_lower: f64,
    pub debug_prev_row_upper: f64,

    reduction_values: HighsDataStack,
    reductions: Vec<(ReductionType, usize)>,
    orig_col_index: Vec<HighsInt>,
    orig_row_index: Vec<HighsInt>,
    /// HiGHS 1.16 (#3218): `origRowType` (1.16.0 line 345), parallel to `orig_row_index`.
    orig_row_type: Vec<OrigRowType>,
    linearly_transformable: Vec<u8>,

    row_values: Vec<Nonzero>,
    col_values: Vec<Nonzero>,
    orig_num_col: HighsInt,
    orig_num_row: HighsInt,
    /// HiGHS 1.16 (#3218): `numAppendedRows`, `nextRowIndex`, `nextColIndex` (1.16.0 lines
    /// 352-354). The 1.15.1 paths keep `next_row_index == orig_num_row` and
    /// `next_col_index == orig_num_col`; only rows/columns appended under
    /// `mip_presolve_row_append` move them apart.
    num_appended_rows: HighsInt,
    next_row_index: HighsInt,
    next_col_index: HighsInt,
}

impl Default for HighsPostsolveStack {
    fn default() -> Self {
        Self {
            debug_prev_numreductions: 0,
            debug_prev_col_lower: 0.0,
            debug_prev_col_upper: 0.0,
            debug_prev_row_lower: 0.0,
            debug_prev_row_upper: 0.0,
            reduction_values: HighsDataStack::default(),
            reductions: Vec::new(),
            orig_col_index: Vec::new(),
            orig_row_index: Vec::new(),
            orig_row_type: Vec::new(),
            linearly_transformable: Vec::new(),
            row_values: Vec::new(),
            col_values: Vec::new(),
            orig_num_col: -1,
            orig_num_row: -1,
            num_appended_rows: 0,
            next_row_index: -1,
            next_col_index: -1,
        }
    }
}

/// `rowValues.emplace_back(origColIndex[rowVal.index()], rowVal.value())`
/// over a slice (shared by the push templates).
fn collect_nonzeros(
    out: &mut Vec<Nonzero>,
    orig_index: &[HighsInt],
    vec: impl IntoIterator<Item = HighsSliceNonzero>,
) {
    out.clear();
    for v in vec {
        out.push(Nonzero::new(orig_index[v.index() as usize], v.value()));
    }
}

impl HighsPostsolveStack {
    /// C5 read accessor (no HiGHS counterpart): the stack's entries tallied by reduction type, in
    /// first-seen order. Reporting only.
    pub fn reduction_counts(&self) -> Vec<(String, i64)> {
        let mut out: Vec<(String, i64)> = Vec::new();
        for (ty, _) in &self.reductions {
            let name = format!("{ty:?}");
            let name = name.strip_prefix('k').unwrap_or(&name).to_string();
            match out.iter_mut().find(|(n, _)| *n == name) {
                Some((_, cnt)) => *cnt += 1,
                None => out.push((name, 1)),
            }
        }
        out
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::reductionAdded
    fn reduction_added(&mut self, type_: ReductionType) {
        let position = self.reduction_values.get_current_data_size();
        self.reductions.push((type_, position));
    }

    /// The recorded reduction types, in order (no HiGHS counterpart: read by
    /// `rusolve presolve-dump`, which diffs presolve against the C++ oracle).
    pub fn debug_reduction_types(&self) -> impl Iterator<Item = ReductionType> + '_ {
        self.reductions.iter().map(|r| r.0)
    }

    /// The indices recorded by reduction `i` for the kinds whose kept/removed choice a
    /// dump diff needs (DuplicateRow, RedundantRow, FixedCol), else an empty string. No
    /// HiGHS counterpart: `rusolve presolve-dump` with `PD_DETAIL` set; the C++ dumper
    /// pops the same records from `reductionValues`.
    pub fn debug_reduction_detail(&self, i: usize) -> String {
        let (type_, position) = self.reductions[i];
        let mut values = self.reduction_values.clone();
        values.set_position(position);
        match type_ {
            ReductionType::kDuplicateRow => {
                let mut d = DuplicateRow::default();
                values.pop(&mut d);
                format!(
                    " row={} dup={} scale={} lt={} ut={}",
                    d.row,
                    d.duplicate_row,
                    d.duplicate_row_scale,
                    d.row_lower_tightened as i32,
                    d.row_upper_tightened as i32
                )
            }
            ReductionType::kRedundantRow => {
                let mut d = RedundantRow::default();
                values.pop(&mut d);
                format!(" row={}", d.row)
            }
            ReductionType::kFixedCol => {
                let mut nz: Vec<Nonzero> = Vec::new();
                values.pop_vec(&mut nz);
                let mut d = FixedCol::default();
                values.pop(&mut d);
                format!(" col={} val={}", d.col, d.fix_value)
            }
            ReductionType::kFreeColSubstitution => {
                let mut nz: Vec<Nonzero> = Vec::new();
                values.pop_vec(&mut nz);
                values.pop_vec(&mut nz);
                let mut d = FreeColSubstitution::default();
                values.pop(&mut d);
                format!(" row={} col={}", d.row, d.col)
            }
            ReductionType::kDoubletonEquation => {
                let mut nz: Vec<Nonzero> = Vec::new();
                values.pop_vec(&mut nz);
                let mut d = DoubletonEquation::default();
                values.pop(&mut d);
                format!(" row={} col={} subst={}", d.row, d.col, d.col_subst)
            }
            ReductionType::kEqualityRowAdditions => {
                let mut nz: Vec<Nonzero> = Vec::new();
                values.pop_vec(&mut nz);
                let targets: Vec<String> = nz.iter().map(|t| t.index.to_string()).collect();
                values.pop_vec(&mut nz);
                let mut d = EqualityRowAdditions::default();
                values.pop(&mut d);
                format!(" row={} targets={}", d.added_eq_row, targets.join(","))
            }
            ReductionType::kForcingRow => {
                let mut nz: Vec<Nonzero> = Vec::new();
                values.pop_vec(&mut nz);
                let mut d = ForcingRow::default();
                values.pop(&mut d);
                format!(" row={}", d.row)
            }
            ReductionType::kForcingColumnRemovedRow => {
                let mut nz: Vec<Nonzero> = Vec::new();
                values.pop_vec(&mut nz);
                let mut d = ForcingColumnRemovedRow::default();
                values.pop(&mut d);
                format!(" row={}", d.row)
            }
            _ => String::new(),
        }
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::getOrigRowsIndex
    pub fn get_orig_rows_index(&self) -> &[HighsInt] {
        &self.orig_row_index
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::getOrigColsIndex
    pub fn get_orig_cols_index(&self) -> &[HighsInt] {
        &self.orig_col_index
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::getOrigRowIndex
    pub fn get_orig_row_index(&self, row: HighsInt) -> HighsInt {
        debug_assert!((row as usize) < self.orig_row_index.len());
        self.orig_row_index[row as usize]
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::getOrigColIndex
    pub fn get_orig_col_index(&self, col: HighsInt) -> HighsInt {
        debug_assert!((col as usize) < self.orig_col_index.len());
        self.orig_col_index[col as usize]
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::appendCutsToModel
    ///
    /// The 1.15.1 text. `orig_row_type` and `next_row_index` are kept in step (inert).
    pub fn append_cuts_to_model(&mut self, num_cuts: HighsInt) {
        let curr_num_row = self.orig_row_index.len();
        let new_num_row = curr_num_row + num_cuts as usize;
        self.orig_row_index.resize(new_num_row, 0);
        self.orig_row_type.resize(new_num_row, OrigRowType::kCut);
        for i in curr_num_row..new_num_row {
            self.orig_row_index[i] = self.orig_num_row;
            self.orig_num_row += 1;
        }
        self.next_row_index = self.orig_num_row;
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::removeCutsFromModel
    pub fn remove_cuts_from_model(&mut self, num_cuts: HighsInt) {
        self.orig_num_row -= num_cuts;

        let mut orig_row_index_size = self.orig_row_index.len();
        let mut i = self.orig_row_index.len();
        while i > 0 {
            if self.orig_row_index[i - 1] < self.orig_num_row {
                break;
            }
            orig_row_index_size -= 1;
            i -= 1;
        }

        self.orig_row_index.resize(orig_row_index_size, 0);
        self.orig_row_type
            .resize(orig_row_index_size, OrigRowType::kOriginal);
        self.next_row_index = self.orig_num_row;
    }

    /// HiGHS 1.16 (#3218): HighsPostsolveStack.h::HighsPostsolveStack::appendToModel
    /// (1.16.0 lines 467-477). C++ passes `HighsInt& numRows` (origNumRow or
    /// numAppendedRows); here `add_to_orig` selects it.
    fn append_to_model(
        &mut self,
        add_to_orig: bool,
        num_rows_to_append: HighsInt,
        row_type: OrigRowType,
    ) {
        if num_rows_to_append <= 0 {
            return;
        }
        let curr_num_row = self.orig_row_index.len();
        let new_num_row = curr_num_row + num_rows_to_append as usize;
        self.orig_row_index.resize(new_num_row, 0);
        self.orig_row_type.resize(new_num_row, row_type);
        for i in curr_num_row..new_num_row {
            self.orig_row_index[i] = self.next_row_index;
            self.next_row_index += 1;
        }
        if add_to_orig {
            self.orig_num_row += num_rows_to_append;
        } else {
            self.num_appended_rows += num_rows_to_append;
        }
    }

    /// HiGHS 1.16 (#3218): HighsPostsolveStack.h::HighsPostsolveStack::appendCutsToModel
    /// (1.16.0 lines 509-511), under `mip_presolve_row_append`.
    pub fn append_cuts_to_model_x1b(&mut self, num_cuts: HighsInt) {
        self.append_to_model(true, num_cuts, OrigRowType::kCut);
    }

    /// HiGHS 1.16 (#3218): HighsPostsolveStack.h::HighsPostsolveStack::appendRowsToModel
    /// (1.16.0 lines 513-515)
    pub fn append_rows_to_model(&mut self, num_rows: HighsInt) {
        self.append_to_model(false, num_rows, OrigRowType::kAppended);
    }

    /// HiGHS 1.16 (#3250): HighsPostsolveStack.h::HighsPostsolveStack::removeCutsFromModel
    /// (1.16.0 lines 517-532), under `mip_presolve_row_append`. `nextRowIndex` is not
    /// decremented.
    pub fn remove_cuts_from_model_x1b(&mut self, num_cuts: HighsInt) {
        if num_cuts <= 0 {
            return;
        }
        self.orig_num_row -= num_cuts;
        let mut new_size = 0usize;
        for i in 0..self.orig_row_index.len() {
            if self.orig_row_type[i] != OrigRowType::kCut {
                if i != new_size {
                    self.orig_row_index[new_size] = self.orig_row_index[i];
                    self.orig_row_type[new_size] = self.orig_row_type[i];
                }
                new_size += 1;
            }
        }
        self.orig_row_index.resize(new_size, 0);
        self.orig_row_type.resize(new_size, OrigRowType::kOriginal);
    }

    /// HiGHS 1.16 (#3218): HighsPostsolveStack.h::HighsPostsolveStack::appendColToModel
    /// (1.16.0 lines 542-545)
    pub fn append_col_to_model(&mut self) {
        self.orig_col_index.push(self.next_col_index);
        self.next_col_index += 1;
        self.linearly_transformable.push(0);
    }

    /// HiGHS 1.16: HighsPostsolveStack.h::HighsPostsolveStack::isOrigCol (line 438)
    pub fn is_orig_col(&self, col: HighsInt) -> bool {
        self.orig_col_index[col as usize] < self.orig_num_col
    }

    /// HiGHS 1.16 (#3250): HighsPostsolveStack.h::HighsPostsolveStack::getOrigCols (lines 441-446)
    pub fn get_orig_cols(&self) -> Vec<HighsInt> {
        (0..self.orig_col_index.len() as HighsInt)
            .filter(|&i| self.is_orig_col(i))
            .collect()
    }

    /// HiGHS 1.16: HighsPostsolveStack.h::HighsPostsolveStack::isOrigRow (lines 448-450)
    pub fn is_orig_row(&self, row: HighsInt) -> bool {
        self.orig_row_type[row as usize] == OrigRowType::kOriginal
    }

    /// HiGHS 1.16: HighsPostsolveStack.h::HighsPostsolveStack::isAppendedRow (lines 452-454)
    pub fn is_appended_row(&self, row: HighsInt) -> bool {
        self.orig_row_type[row as usize] == OrigRowType::kAppended
    }

    /// HiGHS 1.16: HighsPostsolveStack.h::HighsPostsolveStack::isCutRow (lines 456-458)
    pub fn is_cut_row(&self, row: HighsInt) -> bool {
        self.orig_row_type[row as usize] == OrigRowType::kCut
    }

    /// HiGHS 1.16 (#3250): HighsPostsolveStack.h::HighsPostsolveStack::setRowType (line 460)
    pub fn set_row_type(&mut self, row: HighsInt, type_: OrigRowType) {
        self.orig_row_type[row as usize] = type_;
    }

    /// HiGHS 1.16 (#3250): HighsPostsolveStack.h::HighsPostsolveStack::getOrigRows (lines 463-468)
    pub fn get_orig_rows(&self) -> Vec<HighsInt> {
        (0..self.orig_row_type.len() as HighsInt)
            .filter(|&i| self.is_orig_row(i))
            .collect()
    }

    /// HiGHS 1.16 (#3250): HighsPostsolveStack.h::HighsPostsolveStack::getCutRows (lines 471-478),
    /// in reverse order ("keep reverse loop to avoid behavior changes").
    pub fn get_cut_rows(&self) -> Vec<HighsInt> {
        (0..self.orig_row_type.len() as HighsInt)
            .rev()
            .filter(|&i| self.is_cut_row(i))
            .collect()
    }

    /// HiGHS 1.16 (#3250): HighsPostsolveStack.h::HighsPostsolveStack::getNonCutRows (lines 482-487)
    pub fn get_non_cut_rows(&self) -> Vec<HighsInt> {
        (0..self.orig_row_type.len() as HighsInt)
            .filter(|&i| !self.is_cut_row(i))
            .collect()
    }

    /// HiGHS 1.16: HighsPostsolveStack.h::HighsPostsolveStack::getNextRowIndex (line 538)
    pub fn get_next_row_index(&self) -> HighsInt {
        self.next_row_index
    }

    /// HiGHS 1.16: HighsPostsolveStack.h::HighsPostsolveStack::getNextColIndex (line 540)
    pub fn get_next_col_index(&self) -> HighsInt {
        self.next_col_index
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::getOrigNumRow
    pub fn get_orig_num_row(&self) -> HighsInt {
        self.orig_num_row
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::getOrigNumCol
    pub fn get_orig_num_col(&self) -> HighsInt {
        self.orig_num_col
    }

    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::initializeIndexMaps (lines 20-32, 13 lines)
    pub fn initialize_index_maps(&mut self, num_row: HighsInt, num_col: HighsInt) {
        self.orig_num_row = num_row;
        self.orig_num_col = num_col;
        // HiGHS 1.16 (#3218): HighsPostsolveStack.cpp lines 25-26, 31 (1.16.0)
        self.next_row_index = num_row;
        self.next_col_index = num_col;

        self.orig_row_index.resize(num_row as usize, 0);
        for (i, v) in self.orig_row_index.iter_mut().enumerate() {
            *v = i as HighsInt;
        }
        self.orig_row_type
            .resize(num_row as usize, OrigRowType::kOriginal);

        self.orig_col_index.resize(num_col as usize, 0);
        for (i, v) in self.orig_col_index.iter_mut().enumerate() {
            *v = i as HighsInt;
        }

        self.linearly_transformable.resize(num_col as usize, 1);
    }

    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::compressIndexMaps (lines 34-58, 25 lines)
    pub fn compress_index_maps(&mut self, new_row_index: &[HighsInt], new_col_index: &[HighsInt]) {
        let mut num_row = self.orig_row_index.len() as HighsInt;
        for i in 0..new_row_index.len() {
            if new_row_index[i] == -1 {
                num_row -= 1;
            } else {
                self.orig_row_index[new_row_index[i] as usize] = self.orig_row_index[i];
            }
        }
        self.orig_row_index.resize(num_row as usize, 0);
        // HiGHS 1.16 (#3218): compressRowIndexMap also compresses origRowType
        // (HighsPostsolveStack.cpp lines 46-50, 1.16.0), which has the same length.
        let mut num_row_type = self.orig_row_type.len() as HighsInt;
        for i in 0..new_row_index.len() {
            if new_row_index[i] == -1 {
                num_row_type -= 1;
            } else {
                self.orig_row_type[new_row_index[i] as usize] = self.orig_row_type[i];
            }
        }
        self.orig_row_type
            .resize(num_row_type as usize, OrigRowType::kOriginal);

        let mut num_col = self.orig_col_index.len() as HighsInt;
        for i in 0..new_col_index.len() {
            if new_col_index[i] == -1 {
                num_col -= 1;
            } else {
                self.orig_col_index[new_col_index[i] as usize] = self.orig_col_index[i];
            }
        }
        self.orig_col_index.resize(num_col as usize, 0);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::linearTransform
    pub fn linear_transform(&mut self, col: HighsInt, scale: f64, constant: f64) {
        self.reduction_values.push(&LinearTransform {
            scale,
            constant,
            col: self.orig_col_index[col as usize],
        });
        self.reduction_added(ReductionType::kLinearTransform);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::freeColSubstitution
    #[allow(clippy::too_many_arguments)]
    pub fn free_col_substitution(
        &mut self,
        row: HighsInt,
        col: HighsInt,
        rhs: f64,
        col_cost: f64,
        row_type: RowType,
        row_vec: impl IntoIterator<Item = HighsSliceNonzero>,
        col_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        collect_nonzeros(&mut self.row_values, &self.orig_col_index, row_vec);
        collect_nonzeros(&mut self.col_values, &self.orig_row_index, col_vec);

        self.reduction_values.push(&FreeColSubstitution {
            rhs,
            col_cost,
            row: self.orig_row_index[row as usize],
            col: self.orig_col_index[col as usize],
            row_type,
        });
        self.reduction_values.push_vec(&self.row_values);
        self.reduction_values.push_vec(&self.col_values);
        self.reduction_added(ReductionType::kFreeColSubstitution);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::slackColSubstitution
    pub fn slack_col_substitution(
        &mut self,
        row: HighsInt,
        col: HighsInt,
        rhs: f64,
        row_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        collect_nonzeros(&mut self.row_values, &self.orig_col_index, row_vec);

        self.reduction_values.push(&SlackColSubstitution {
            rhs,
            row: self.orig_row_index[row as usize],
            col: self.orig_col_index[col as usize],
        });
        self.reduction_values.push_vec(&self.row_values);
        self.reduction_added(ReductionType::kSlackColSubstitution);
    }

    /// HiGHS 1.16 (#2962): HighsPostsolveStack.h::HighsPostsolveStack::zeroCostSingleton (lines 597-610)
    #[allow(clippy::too_many_arguments)]
    pub fn zero_cost_singleton(
        &mut self,
        row: HighsInt,
        col: HighsInt,
        orig_row_lower: f64,
        orig_row_upper: f64,
        new_row_lb: f64,
        new_row_ub: f64,
        lb: f64,
        ub: f64,
        coef: f64,
        row_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        collect_nonzeros(&mut self.row_values, &self.orig_col_index, row_vec);

        self.reduction_values.push(&ZeroObjSingletonContinuousCol {
            orig_row_lower,
            orig_row_upper,
            new_row_lb,
            new_row_ub,
            lb,
            ub,
            coef,
            col: self.orig_col_index[col as usize],
            row: self.orig_row_index[row as usize],
        });
        self.reduction_values.push_vec(&self.row_values);
        self.reduction_added(ReductionType::kZeroObjSingletonContinuousCol);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::doubletonEquation
    #[allow(clippy::too_many_arguments)]
    pub fn doubleton_equation(
        &mut self,
        row: HighsInt,
        col_subst: HighsInt,
        col: HighsInt,
        coef_subst: f64,
        coef: f64,
        rhs: f64,
        subst_lower: f64,
        subst_upper: f64,
        subst_cost: f64,
        lower_tightened: bool,
        upper_tightened: bool,
        row_type: RowType,
        col_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        collect_nonzeros(&mut self.col_values, &self.orig_row_index, col_vec);

        self.reduction_values.push(&DoubletonEquation {
            coef,
            coef_subst,
            rhs,
            subst_lower,
            subst_upper,
            subst_cost,
            row: if row == -1 {
                -1
            } else {
                self.orig_row_index[row as usize]
            },
            col_subst: self.orig_col_index[col_subst as usize],
            col: self.orig_col_index[col as usize],
            lower_tightened,
            upper_tightened,
            row_type,
        });
        self.reduction_values.push_vec(&self.col_values);
        self.reduction_added(ReductionType::kDoubletonEquation);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::equalityRowAddition
    pub fn equality_row_addition(
        &mut self,
        row: HighsInt,
        added_eq_row: HighsInt,
        eq_row_scale: f64,
        eq_row_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        collect_nonzeros(&mut self.row_values, &self.orig_col_index, eq_row_vec);

        self.reduction_values.push(&EqualityRowAddition {
            row: self.orig_row_index[row as usize],
            added_eq_row: self.orig_row_index[added_eq_row as usize],
            eq_row_scale,
        });
        self.reduction_values.push_vec(&self.row_values);
        self.reduction_added(ReductionType::kEqualityRowAddition);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::equalityRowAdditions
    /// HiGHS 1.16 (#3324): HighsPostsolveStack.h::equalityRowAdditions (lines 763-782)
    /// `h16_eq_row_adds` is passed in because the postsolve stack has no options in
    /// reach at push time; `false` reproduces 1.15.1.
    pub fn equality_row_additions(
        &mut self,
        added_eq_row: HighsInt,
        eq_row_vec: impl IntoIterator<Item = HighsSliceNonzero>,
        target_rows: &[Nonzero],
        h16_eq_row_adds: bool,
    ) {
        collect_nonzeros(&mut self.row_values, &self.orig_col_index, eq_row_vec);

        self.reduction_values.push(&EqualityRowAdditions {
            added_eq_row: self.orig_row_index[added_eq_row as usize],
        });
        self.reduction_values.push_vec(&self.row_values);
        // HiGHS 1.16 #3324 (h16_eq_row_adds): HighsPostsolveStack.h lines 772-779 (1.16.0)
        if h16_eq_row_adds {
            // targetRows are in the presolved space; postsolve needs the original space.
            let orig_target_rows: Vec<Nonzero> = target_rows
                .iter()
                .map(|t| Nonzero::new(self.orig_row_index[t.index as usize], t.value))
                .collect();
            self.reduction_values.push_vec(&orig_target_rows);
        } else {
            self.reduction_values.push_vec(target_rows);
        }
        self.reduction_added(ReductionType::kEqualityRowAdditions);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::singletonRow
    pub fn singleton_row(
        &mut self,
        row: HighsInt,
        col: HighsInt,
        coef: f64,
        tightened_col_lower: bool,
        tightened_col_upper: bool,
    ) {
        self.reduction_values.push(&SingletonRow {
            coef,
            row: self.orig_row_index[row as usize],
            col: self.orig_col_index[col as usize],
            col_lower_tightened: tightened_col_lower,
            col_upper_tightened: tightened_col_upper,
        });
        self.reduction_added(ReductionType::kSingletonRow);
    }

    fn push_fixed_col(
        &mut self,
        col: HighsInt,
        fix_value: f64,
        col_cost: f64,
        fix_type: HighsBasisStatus,
        col_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        collect_nonzeros(&mut self.col_values, &self.orig_row_index, col_vec);

        self.reduction_values.push(&FixedCol {
            fix_value,
            col_cost,
            col: self.orig_col_index[col as usize],
            fix_type,
        });
        self.reduction_values.push_vec(&self.col_values);
        self.reduction_added(ReductionType::kFixedCol);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::fixedColAtLower
    pub fn fixed_col_at_lower(
        &mut self,
        col: HighsInt,
        fix_value: f64,
        col_cost: f64,
        col_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        debug_assert!(fix_value.is_finite());
        self.push_fixed_col(col, fix_value, col_cost, HighsBasisStatus::kLower, col_vec);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::fixedColAtUpper
    pub fn fixed_col_at_upper(
        &mut self,
        col: HighsInt,
        fix_value: f64,
        col_cost: f64,
        col_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        debug_assert!(fix_value.is_finite());
        self.push_fixed_col(col, fix_value, col_cost, HighsBasisStatus::kUpper, col_vec);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::fixedColAtZero
    pub fn fixed_col_at_zero(
        &mut self,
        col: HighsInt,
        col_cost: f64,
        col_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        self.push_fixed_col(col, 0.0, col_cost, HighsBasisStatus::kZero, col_vec);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::removedFixedCol
    pub fn removed_fixed_col(
        &mut self,
        col: HighsInt,
        fix_value: f64,
        col_cost: f64,
        col_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        debug_assert!(fix_value.is_finite());
        self.push_fixed_col(
            col,
            fix_value,
            col_cost,
            HighsBasisStatus::kNonbasic,
            col_vec,
        );
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::redundantRow
    pub fn redundant_row(&mut self, row: HighsInt) {
        self.reduction_values.push(&RedundantRow {
            row: self.orig_row_index[row as usize],
        });
        self.reduction_added(ReductionType::kRedundantRow);
    }

    /// HiGHS 1.16 (#3218): HighsPostsolveStack.h::HighsPostsolveStack::impliedEquation (lines 864-874)
    pub fn implied_equation(
        &mut self,
        row: HighsInt,
        at_lower: bool,
        row_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        collect_nonzeros(&mut self.row_values, &self.orig_col_index, row_vec);
        self.reduction_values.push(&ImpliedEquation {
            row: self.orig_row_index[row as usize],
            at_lower,
        });
        self.reduction_values.push_vec(&self.row_values);
        self.reduction_added(ReductionType::kImpliedEquation);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::forcingRow
    pub fn forcing_row(
        &mut self,
        row: HighsInt,
        row_vec: impl IntoIterator<Item = HighsSliceNonzero>,
        side: f64,
        row_type: RowType,
    ) {
        collect_nonzeros(&mut self.row_values, &self.orig_col_index, row_vec);

        self.reduction_values.push(&ForcingRow {
            side,
            row: self.orig_row_index[row as usize],
            row_type,
        });
        self.reduction_values.push_vec(&self.row_values);
        self.reduction_added(ReductionType::kForcingRow);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::forcingColumn
    pub fn forcing_column(
        &mut self,
        col: HighsInt,
        col_vec: impl IntoIterator<Item = HighsSliceNonzero>,
        cost: f64,
        bound_val: f64,
        at_infinite_upper: bool,
        col_integral: bool,
    ) {
        collect_nonzeros(&mut self.col_values, &self.orig_row_index, col_vec);

        self.reduction_values.push(&ForcingColumn {
            col_cost: cost,
            col_bound: bound_val,
            col: self.orig_col_index[col as usize],
            at_infinite_upper,
            col_integral,
        });
        self.reduction_values.push_vec(&self.col_values);
        self.reduction_added(ReductionType::kForcingColumn);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::forcingColumnRemovedRow
    pub fn forcing_column_removed_row(
        &mut self,
        forcing_col: HighsInt,
        row: HighsInt,
        rhs: f64,
        row_vec: impl IntoIterator<Item = HighsSliceNonzero>,
    ) {
        self.row_values.clear();
        for row_val in row_vec {
            if row_val.index() != forcing_col {
                self.row_values.push(Nonzero::new(
                    self.orig_col_index[row_val.index() as usize],
                    row_val.value(),
                ));
            }
        }

        self.reduction_values.push(&ForcingColumnRemovedRow {
            rhs,
            row: self.orig_row_index[row as usize],
        });
        self.reduction_values.push_vec(&self.row_values);
        self.reduction_added(ReductionType::kForcingColumnRemovedRow);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::duplicateRow
    pub fn duplicate_row(
        &mut self,
        row: HighsInt,
        row_upper_tightened: bool,
        row_lower_tightened: bool,
        duplicate_row: HighsInt,
        duplicate_row_scale: f64,
    ) {
        self.reduction_values.push(&DuplicateRow {
            duplicate_row_scale,
            duplicate_row: self.orig_row_index[duplicate_row as usize],
            row: self.orig_row_index[row as usize],
            row_lower_tightened,
            row_upper_tightened,
        });
        self.reduction_added(ReductionType::kDuplicateRow);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::duplicateColumn
    #[allow(clippy::too_many_arguments)]
    pub fn duplicate_column(
        &mut self,
        col_scale: f64,
        col_lower: f64,
        col_upper: f64,
        duplicate_col_lower: f64,
        duplicate_col_upper: f64,
        col: HighsInt,
        duplicate_col: HighsInt,
        col_integral: bool,
        duplicate_col_integral: bool,
        ok_merge_tolerance: f64,
    ) -> bool {
        let orig_col = self.orig_col_index[col as usize];
        let orig_duplicate_col = self.orig_col_index[duplicate_col as usize];
        let debug_values = DuplicateColumn {
            col_scale,
            col_lower,
            col_upper,
            duplicate_col_lower,
            duplicate_col_upper,
            col: orig_col,
            duplicate_col: orig_duplicate_col,
            col_integral,
            duplicate_col_integral,
        };
        let ok_merge = debug_values.ok_merge(ok_merge_tolerance);
        let prevent_illegal_merge = true;
        if !ok_merge && prevent_illegal_merge {
            return false;
        }
        self.reduction_values.push(&debug_values);

        self.reduction_added(ReductionType::kDuplicateColumn);

        // mark columns as not linearly transformable
        self.linearly_transformable[orig_col as usize] = 0;
        self.linearly_transformable[orig_duplicate_col as usize] = 0;
        true
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::getReducedPrimalSolution (lines 541-573, 33 lines)
    pub fn get_reduced_primal_solution(&mut self, orig_primal_solution: &[f64]) -> Vec<f64> {
        let mut reduced_solution = orig_primal_solution.to_vec();
        // HiGHS 1.16 (#3218): `reducedSolution.resize(nextColIndex, 0.0)` (1.16.0 line 962).
        // Columns are appended only under `mip_presolve_row_append` (FM objective
        // reformulation); otherwise next_col_index == orig_num_col and this is skipped, as in
        // 1.15.1. 1.16.0 leaves the appended column at 0 (no
        // FourierMotzkinObjCol::transformToPresolvedSpace call).
        if self.next_col_index != self.orig_num_col {
            reduced_solution.resize(self.next_col_index as usize, 0.0);
        }

        for primal_col_transformation in self.reductions.iter() {
            match primal_col_transformation.0 {
                ReductionType::kDuplicateColumn => {
                    let mut duplicate_col_reduction = DuplicateColumn::default();
                    self.reduction_values
                        .set_position(primal_col_transformation.1);
                    self.reduction_values.pop(&mut duplicate_col_reduction);
                    duplicate_col_reduction.transform_to_presolved_space(&mut reduced_solution);
                }
                ReductionType::kLinearTransform => {
                    self.reduction_values
                        .set_position(primal_col_transformation.1);
                    let mut linear_transform = LinearTransform::default();
                    self.reduction_values.pop(&mut linear_transform);
                    linear_transform.transform_to_presolved_space(&mut reduced_solution);
                }
                _ => continue,
            }
        }

        let reduced_num_col = self.orig_col_index.len();
        for i in 0..reduced_num_col {
            let value = reduced_solution[self.orig_col_index[i] as usize];
            reduced_solution[i] = value;
        }

        reduced_solution.resize(reduced_num_col, 0.0);
        reduced_solution
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::isColLinearlyTransformable
    pub fn is_col_linearly_transformable(&self, col: HighsInt) -> bool {
        debug_assert!(col >= 0);
        debug_assert!((col as usize) < self.orig_col_index.len());
        self.linearly_transformable[self.orig_col_index[col as usize] as usize] != 0
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::undoIterateBackwards (lines 581-601, 21 lines)
    ///
    /// Port the `#else` (non-`DEBUG_EXTRA`) branch. `values.resize(origSize)`
    /// -> `values.resize(orig_size as usize, T::default())`.
    pub fn undo_iterate_backwards<T: Clone + Default>(
        values: &mut Vec<T>,
        index: &[HighsInt],
        orig_size: HighsInt,
    ) {
        values.resize(orig_size as usize, T::default());
        let mut i = index.len();
        while i > 0 {
            debug_assert!((index[i - 1] as usize) >= i - 1);
            let value = values[i - 1].clone();
            values[index[i - 1] as usize] = value;
            i -= 1;
        }
    }

    /// HiGHS 1.16 (#3218): HighsPostsolveStack.h::HighsPostsolveStack::undoIterateBackwards
    /// (1.16.0 lines 1000-1025, `#else` branch): fills the new tail and every vacated slot
    /// with `zero`. Used by `undo` under `mip_presolve_row_append`.
    pub fn undo_iterate_backwards_h16<T: Clone + Default>(
        values: &mut Vec<T>,
        index: &[HighsInt],
        orig_size: HighsInt,
        zero: T,
    ) {
        values.resize(orig_size as usize, T::default());
        for v in values.iter_mut().take(orig_size as usize).skip(index.len()) {
            *v = zero.clone();
        }
        let mut i = index.len();
        while i > 0 {
            let to_i = index[i - 1] as usize;
            debug_assert!(to_i >= i - 1);
            debug_assert!(to_i < orig_size as usize);
            let value = values[i - 1].clone();
            values[to_i] = value;
            if to_i > i - 1 {
                values[i - 1] = zero.clone();
            }
            i -= 1;
        }
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::containsNanOrInf
    pub fn contains_nan_or_inf(&self, v: &[f64]) -> bool {
        v.iter().any(|d| d.is_nan() || d.is_infinite())
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::undo (lines 611-777, 167 lines)
    ///
    /// C++ default arguments: `report_col = -1, thread_safe = false`.
    /// See the module doc for the `thread_safe` reference binding.
    pub fn undo(
        &mut self,
        options: &HighsOptions,
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
        report_col: HighsInt,
        thread_safe: bool,
    ) {
        // HighsDataStack reductionValuesCopy; std::vector<Nonzero>
        // colValuesCopy, rowValuesCopy; and the references
        // reductionValues_/colValues_/rowValues_ bound to the copies or the
        // members (see the module doc): taken out of self when !thread_safe
        // and moved back before returning.
        let (mut reduction_values_, mut col_values_, mut row_values_) = if thread_safe {
            let mut reduction_values_copy = self.reduction_values.clone();
            reduction_values_copy.reset_position();
            (
                reduction_values_copy,
                self.col_values.clone(),
                self.row_values.clone(),
            )
        } else {
            self.reduction_values.reset_position();
            (
                std::mem::take(&mut self.reduction_values),
                std::mem::take(&mut self.col_values),
                std::mem::take(&mut self.row_values),
            )
        };

        // Verify that undo can be performed
        debug_assert!(solution.value_valid);
        let perform_dual_postsolve = solution.dual_valid;
        let perform_basis_postsolve = basis.valid;

        // HiGHS #3218 (mip_presolve_row_append): HighsPostsolveStack.h::undo (1.16.0 lines
        // 1058-1081): expand to nextColIndex / nextRowIndex with zero fill, and resize back
        // to origNumCol / origNumRow after the loop (lines 1300-1306). The isModelRow guards
        // of the undo functions (removed by #3218) then always hold, since every recorded
        // row index is below nextRowIndex.
        let x1b_row_append = options.mip_presolve_row_append;
        if x1b_row_append {
            debug_assert!(self.next_col_index > 0);
            let (nc, nr) = (self.next_col_index, self.next_row_index);
            Self::undo_iterate_backwards_h16(
                &mut solution.col_value,
                &self.orig_col_index,
                nc,
                0.0,
            );
            debug_assert!(nr >= 0);
            Self::undo_iterate_backwards_h16(
                &mut solution.row_value,
                &self.orig_row_index,
                nr,
                0.0,
            );
            if perform_dual_postsolve {
                Self::undo_iterate_backwards_h16(
                    &mut solution.col_dual,
                    &self.orig_col_index,
                    nc,
                    0.0,
                );
                Self::undo_iterate_backwards_h16(
                    &mut solution.row_dual,
                    &self.orig_row_index,
                    nr,
                    0.0,
                );
            }
            if perform_basis_postsolve {
                Self::undo_iterate_backwards_h16(
                    &mut basis.col_status,
                    &self.orig_col_index,
                    nc,
                    HighsBasisStatus::kNonbasic,
                );
                Self::undo_iterate_backwards_h16(
                    &mut basis.row_status,
                    &self.orig_row_index,
                    nr,
                    HighsBasisStatus::kNonbasic,
                );
            }
        }

        // expand solution to original index space
        debug_assert!(self.orig_num_col > 0);
        if !x1b_row_append {
            Self::undo_iterate_backwards(
                &mut solution.col_value,
                &self.orig_col_index,
                self.orig_num_col,
            );

            debug_assert!(self.orig_num_row >= 0);
            Self::undo_iterate_backwards(
                &mut solution.row_value,
                &self.orig_row_index,
                self.orig_num_row,
            );

            if perform_dual_postsolve {
                // if dual solution is given, expand dual solution and basis to original
                // index space
                Self::undo_iterate_backwards(
                    &mut solution.col_dual,
                    &self.orig_col_index,
                    self.orig_num_col,
                );

                Self::undo_iterate_backwards(
                    &mut solution.row_dual,
                    &self.orig_row_index,
                    self.orig_num_row,
                );
            }

            if perform_basis_postsolve {
                // if basis is given, expand basis status values to original index space
                Self::undo_iterate_backwards(
                    &mut basis.col_status,
                    &self.orig_col_index,
                    self.orig_num_col,
                );

                Self::undo_iterate_backwards(
                    &mut basis.row_status,
                    &self.orig_row_index,
                    self.orig_num_row,
                );
            }
        }

        // now undo the changes
        let mut i = self.reductions.len();
        while i > 0 {
            if report_col >= 0 {
                println!(
                    "Before  reduction {:2} (type {:2}): col_value[{:2}] = {}",
                    i - 1,
                    self.reductions[i - 1].0 as i32,
                    report_col,
                    solution.col_value[report_col as usize]
                );
            }
            match self.reductions[i - 1].0 {
                ReductionType::kLinearTransform => {
                    let mut reduction = LinearTransform::default();
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, solution);
                }
                ReductionType::kFreeColSubstitution => {
                    let mut reduction = FreeColSubstitution::default();
                    reduction_values_.pop_vec(&mut col_values_);
                    reduction_values_.pop_vec(&mut row_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &row_values_, &col_values_, solution, basis);
                }
                ReductionType::kDoubletonEquation => {
                    let mut reduction = DoubletonEquation::default();
                    reduction_values_.pop_vec(&mut col_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &col_values_, solution, basis);
                }
                ReductionType::kEqualityRowAddition => {
                    let mut reduction = EqualityRowAddition::default();
                    reduction_values_.pop_vec(&mut row_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &row_values_, solution, basis);
                }
                ReductionType::kEqualityRowAdditions => {
                    let mut reduction = EqualityRowAdditions::default();
                    reduction_values_.pop_vec(&mut col_values_);
                    reduction_values_.pop_vec(&mut row_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &row_values_, &col_values_, solution, basis);
                }
                ReductionType::kSingletonRow => {
                    let mut reduction = SingletonRow::default();
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, solution, basis);
                }
                ReductionType::kFixedCol => {
                    let mut reduction = FixedCol::default();
                    reduction_values_.pop_vec(&mut col_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &col_values_, solution, basis);
                }
                ReductionType::kRedundantRow => {
                    let mut reduction = RedundantRow::default();
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, solution, basis);
                }
                ReductionType::kForcingRow => {
                    let mut reduction = ForcingRow::default();
                    reduction_values_.pop_vec(&mut row_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &row_values_, solution, basis);
                }
                ReductionType::kForcingColumn => {
                    let mut reduction = ForcingColumn::default();
                    reduction_values_.pop_vec(&mut col_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &col_values_, solution, basis);
                }
                ReductionType::kForcingColumnRemovedRow => {
                    let mut reduction = ForcingColumnRemovedRow::default();
                    reduction_values_.pop_vec(&mut row_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &row_values_, solution, basis);
                }
                ReductionType::kDuplicateRow => {
                    let mut reduction = DuplicateRow::default();
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, solution, basis);
                }
                ReductionType::kDuplicateColumn => {
                    let mut reduction = DuplicateColumn::default();
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, solution, basis);
                }
                ReductionType::kSlackColSubstitution => {
                    let mut reduction = SlackColSubstitution::default();
                    reduction_values_.pop_vec(&mut row_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &row_values_, solution, basis);
                }
                ReductionType::kZeroObjSingletonContinuousCol => {
                    let mut reduction = ZeroObjSingletonContinuousCol::default();
                    reduction_values_.pop_vec(&mut row_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(options, &row_values_, solution, basis);
                }
                ReductionType::kImpliedEquation => {
                    let mut reduction = ImpliedEquation::default();
                    reduction_values_.pop_vec(&mut row_values_);
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(&row_values_, solution);
                }
                ReductionType::kFourierMotzkinBlock => {
                    let steps = Self::pop_fourier_motzkin_block(&mut reduction_values_);
                    Self::undo_fourier_motzkin_block(&steps, options, solution, basis);
                }
                ReductionType::kFourierMotzkinObjCol => {
                    let mut cost_entries: Vec<Nonzero> = Vec::new();
                    reduction_values_.pop_vec(&mut cost_entries);
                    let mut reduction = FourierMotzkinObjCol::default();
                    reduction_values_.pop(&mut reduction);
                    reduction.undo(&cost_entries, solution);
                }
                // C++ `default:` (the Rust enum is exhaustive, so this arm
                // is statically unreachable).
                #[allow(unreachable_patterns)]
                _ => {
                    println!(
                        "Reduction case {} not handled",
                        self.reductions[i - 1].0 as i32
                    );
                    if kAllowDeveloperAssert {
                        if cfg!(debug_assertions) {
                            panic!("assertion failed: 1 == 0");
                        }
                    }
                }
            }
            i -= 1;
        }
        if report_col >= 0 {
            println!(
                "After last reduction: col_value[{:2}] = {}",
                report_col, solution.col_value[report_col as usize]
            );
        }

        if x1b_row_append {
            // HighsPostsolveStack.h lines 1300-1306 (1.16.0)
            let (nc, nr) = (self.orig_num_col as usize, self.orig_num_row as usize);
            solution.col_value.resize(nc, 0.0);
            if perform_dual_postsolve {
                solution.col_dual.resize(nc, 0.0);
            }
            if perform_basis_postsolve {
                basis.col_status.resize(nc, HighsBasisStatus::default());
            }
            solution.row_value.resize(nr, 0.0);
            if perform_dual_postsolve {
                solution.row_dual.resize(nr, 0.0);
            }
            if perform_basis_postsolve {
                basis.row_status.resize(nr, HighsBasisStatus::default());
            }
        }

        // #ifdef DEBUG_EXTRA (not compiled in the reference build):
        // solution should not contain NaN or Inf
        // assert(!containsNanOrInf(solution.col_value));
        // row values are not determined by postsolve
        // assert(!containsNanOrInf(solution.row_value));
        // assert(!containsNanOrInf(solution.col_dual));
        // assert(!containsNanOrInf(solution.row_dual));

        if !thread_safe {
            self.reduction_values = reduction_values_;
            self.col_values = col_values_;
            self.row_values = row_values_;
        }
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::undoPrimal (lines 780-790, 11 lines)
    ///
    /// C++ default arguments: `report_col = -1, thread_safe = false`.
    /// `HighsBasis basis; basis.valid = false;` -> `HighsBasis::default()`
    /// then set `valid = false` explicitly.
    pub fn undo_primal(
        &mut self,
        options: &HighsOptions,
        solution: &mut HighsSolution,
        report_col: HighsInt,
        thread_safe: bool,
    ) {
        // Call to reductionValues.resetPosition(); seems unnecessary as
        // it's the first thing done in undo
        if !thread_safe {
            self.reduction_values.reset_position();
        }
        let mut basis = HighsBasis::default();
        basis.valid = false;
        solution.dual_valid = false;
        self.undo(options, solution, &mut basis, report_col, thread_safe);
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::undoUntil (lines 806-953, 148 lines)
    ///
    /// The `reductionValues.pop(colValues)` + `reduction.undo(..., colValues, ...)`
    /// pattern borrows two members: take `self.col_values` / `self.row_values`
    /// out with `std::mem::take` around the loop and put them back after it.
    pub fn undo_until(
        &mut self,
        options: &HighsOptions,
        solution: &mut HighsSolution,
        basis: &mut HighsBasis,
        num_reductions: usize,
    ) {
        self.reduction_values.reset_position();

        // Do these returns ever happen? How is it known that undo has not
        // been performed?
        debug_assert!(solution.col_value.len() == self.orig_col_index.len());
        debug_assert!(solution.row_value.len() == self.orig_row_index.len());
        // This should be a better measure of whether undo can be
        // performed
        debug_assert!(solution.value_valid);
        if solution.col_value.len() != self.orig_col_index.len() {
            return;
        }
        if solution.row_value.len() != self.orig_row_index.len() {
            return;
        }

        let perform_dual_postsolve = solution.dual_valid;
        debug_assert!(
            (solution.col_dual.len() == solution.col_value.len()) == perform_dual_postsolve
        );
        let perform_basis_postsolve = basis.valid;

        // expand solution to original index space
        Self::undo_iterate_backwards(
            &mut solution.col_value,
            &self.orig_col_index,
            self.orig_num_col,
        );

        Self::undo_iterate_backwards(
            &mut solution.row_value,
            &self.orig_row_index,
            self.orig_num_row,
        );

        if perform_dual_postsolve {
            // if dual solution is given, expand dual solution and basis to original
            // index space
            Self::undo_iterate_backwards(
                &mut solution.col_dual,
                &self.orig_col_index,
                self.orig_num_col,
            );

            Self::undo_iterate_backwards(
                &mut solution.row_dual,
                &self.orig_row_index,
                self.orig_num_row,
            );
        }

        if perform_basis_postsolve {
            // if basis is given, expand basis status values to original index space
            Self::undo_iterate_backwards(
                &mut basis.col_status,
                &self.orig_col_index,
                self.orig_num_col,
            );

            Self::undo_iterate_backwards(
                &mut basis.row_status,
                &self.orig_row_index,
                self.orig_num_row,
            );
        }

        // now undo the changes. reductionValues stays in self; the two
        // member vectors the pops write into are taken out around the loop
        // (module doc) and moved back after it.
        let mut col_values = std::mem::take(&mut self.col_values);
        let mut row_values = std::mem::take(&mut self.row_values);
        let mut i = self.reductions.len();
        while i > num_reductions {
            match self.reductions[i - 1].0 {
                ReductionType::kLinearTransform => {
                    let mut reduction = LinearTransform::default();
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, solution);
                }
                ReductionType::kFreeColSubstitution => {
                    let mut reduction = FreeColSubstitution::default();
                    self.reduction_values.pop_vec(&mut col_values);
                    self.reduction_values.pop_vec(&mut row_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &row_values, &col_values, solution, basis);
                }
                ReductionType::kDoubletonEquation => {
                    let mut reduction = DoubletonEquation::default();
                    self.reduction_values.pop_vec(&mut col_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &col_values, solution, basis);
                }
                ReductionType::kEqualityRowAddition => {
                    let mut reduction = EqualityRowAddition::default();
                    self.reduction_values.pop_vec(&mut row_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &row_values, solution, basis);
                }
                ReductionType::kEqualityRowAdditions => {
                    let mut reduction = EqualityRowAdditions::default();
                    self.reduction_values.pop_vec(&mut col_values);
                    self.reduction_values.pop_vec(&mut row_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &row_values, &col_values, solution, basis);
                }
                ReductionType::kSingletonRow => {
                    let mut reduction = SingletonRow::default();
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, solution, basis);
                }
                ReductionType::kFixedCol => {
                    let mut reduction = FixedCol::default();
                    self.reduction_values.pop_vec(&mut col_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &col_values, solution, basis);
                }
                ReductionType::kRedundantRow => {
                    let mut reduction = RedundantRow::default();
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, solution, basis);
                }
                ReductionType::kForcingRow => {
                    let mut reduction = ForcingRow::default();
                    self.reduction_values.pop_vec(&mut row_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &row_values, solution, basis);
                }
                ReductionType::kForcingColumn => {
                    let mut reduction = ForcingColumn::default();
                    self.reduction_values.pop_vec(&mut col_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &col_values, solution, basis);
                }
                ReductionType::kForcingColumnRemovedRow => {
                    let mut reduction = ForcingColumnRemovedRow::default();
                    self.reduction_values.pop_vec(&mut row_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &row_values, solution, basis);
                }
                ReductionType::kDuplicateRow => {
                    let mut reduction = DuplicateRow::default();
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, solution, basis);
                }
                ReductionType::kDuplicateColumn => {
                    let mut reduction = DuplicateColumn::default();
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, solution, basis);
                }
                ReductionType::kSlackColSubstitution => {
                    let mut reduction = SlackColSubstitution::default();
                    self.reduction_values.pop_vec(&mut row_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &row_values, solution, basis);
                }
                ReductionType::kZeroObjSingletonContinuousCol => {
                    let mut reduction = ZeroObjSingletonContinuousCol::default();
                    self.reduction_values.pop_vec(&mut row_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(options, &row_values, solution, basis);
                }
                ReductionType::kImpliedEquation => {
                    let mut reduction = ImpliedEquation::default();
                    self.reduction_values.pop_vec(&mut row_values);
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(&row_values, solution);
                }
                ReductionType::kFourierMotzkinBlock => {
                    let steps = Self::pop_fourier_motzkin_block(&mut self.reduction_values);
                    Self::undo_fourier_motzkin_block(&steps, options, solution, basis);
                }
                ReductionType::kFourierMotzkinObjCol => {
                    let mut cost_entries: Vec<Nonzero> = Vec::new();
                    self.reduction_values.pop_vec(&mut cost_entries);
                    let mut reduction = FourierMotzkinObjCol::default();
                    self.reduction_values.pop(&mut reduction);
                    reduction.undo(&cost_entries, solution);
                }
            }
            i -= 1;
        }
        self.col_values = col_values;
        self.row_values = row_values;

        // #ifdef DEBUG_EXTRA (not compiled in the reference build):
        // solution should not contain NaN or Inf
        // assert(!containsNanOrInf(solution.col_value));
        // row values are not determined by postsolve
        // assert(!containsNanOrInf(solution.row_value));
        // assert(!containsNanOrInf(solution.col_dual));
        // assert(!containsNanOrInf(solution.row_dual));
    }

    /// HiGHS: HighsPostsolveStack.h::HighsPostsolveStack::numReductions
    pub fn num_reductions(&self) -> usize {
        self.reductions.len()
    }
}

impl LinearTransform {
    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::LinearTransform::undo (lines 60-66, 7 lines)
    pub fn undo(&self, options: &HighsOptions, solution: &mut HighsSolution) {
        let col = self.col as usize;
        solution.col_value[col] *= self.scale;
        solution.col_value[col] += self.constant;

        if solution.dual_valid {
            solution.col_dual[col] /= self.scale;
        }
    }

    /// HiGHS: HighsPostsolveStack.cpp::HighsPostsolveStack::LinearTransform::transformToPresolvedSpace (lines 68-72, 5 lines)
    pub fn transform_to_presolved_space(&self, primal_sol: &mut [f64]) {
        let col = self.col as usize;
        primal_sol[col] -= self.constant;
        primal_sol[col] /= self.scale;
    }
}

#[cfg(test)]
mod h16_eq_row_adds_tests {
    use super::*;
    use crate::lp_data::highs_options::HighsOptions;

    /// Stack for a 3-row, 2-col model whose row 0 was removed, so presolved row 0 is
    /// original row 1 and presolved row 1 is original row 2. The equality row is presolved
    /// row 0; the target row is presolved row 1 with scale 2.0. Presolved duals are [1, 5].
    fn undo_with_target(h16_eq_row_adds: bool) -> HighsSolution {
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(3, 2);
        stack.compress_index_maps(&[-1, 0, 1], &[0, 1]);
        stack.equality_row_additions(
            0,
            std::iter::empty::<HighsSliceNonzero>(),
            &[Nonzero::new(1, 2.0)],
            h16_eq_row_adds,
        );

        let mut solution = HighsSolution::default();
        solution.value_valid = true;
        solution.dual_valid = true;
        solution.col_value = vec![0.0, 0.0];
        solution.col_dual = vec![0.0, 0.0];
        solution.row_value = vec![0.0, 0.0];
        solution.row_dual = vec![1.0, 5.0];
        let mut basis = HighsBasis::default();
        let options = HighsOptions::default();
        stack.undo(&options, &mut solution, &mut basis, -1, false);
        solution
    }

    #[test]
    fn h16_eq_row_adds_on_maps_target_rows_to_orig_space() {
        // Target row in original space is row 2 (dual 5): row_dual[1] = 1 + 2 * 5 = 11.
        let solution = undo_with_target(true);
        assert_eq!(solution.row_dual[1], 11.0);
    }

    #[test]
    fn h16_eq_row_adds_off_keeps_presolved_target_rows() {
        // 1.15.1 stores the presolved index 1, which is original row 1 (the equality
        // row itself): row_dual[1] = 1 + 2 * 1 = 3.
        let solution = undo_with_target(false);
        assert_eq!(solution.row_dual[1], 3.0);
    }

    /// TestPresolve.cpp "test-sparsify-primal-dual-postsolve" (3324.diff lines 5-97):
    /// the empty first row is removed before sparsify, so presolved and original row
    /// indices differ. Solved with IPX and no crossover; the C++ test only requires kOk
    /// and kOptimal. Sparsify does call equality_row_additions here, but in this port the
    /// row index maps are still the identity at that point, so this test does not
    /// separate on from off: the stack-level tests above do.
    fn sparsify_primal_dual_postsolve(h16_eq_row_adds: bool) -> (HighsStatus, HighsModelStatus) {
        use crate::prelude::*;
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 3;
        lp.col_cost_ = vec![0.0, 0.0, 1.0];
        lp.col_lower_ = vec![-kHighsInf, -kHighsInf, -kHighsInf];
        lp.col_upper_ = vec![1.0, 1.0, 1.0];
        lp.a_matrix_.start_ = vec![0, 2, 4, 6];
        lp.a_matrix_.index_ = vec![1, 2, 1, 2, 1, 2];
        lp.a_matrix_.value_ = vec![8.0, -1.0, 1.0, -1.0, 1.0, -1.0];
        lp.row_lower_ = vec![0.0, -kHighsInf, 0.0];
        lp.row_upper_ = vec![0.0, 0.0, 0.0];
        let mut h = Highs::default();
        h.set_option_value_bool("output_flag", false);
        h.set_option_value_bool("h16_eq_row_adds", h16_eq_row_adds);
        assert_eq!(h.pass_model(lp), HighsStatus::kOk);
        h.set_option_value_str("solver", "ipx");
        h.set_option_value_str("run_crossover", "off");
        let status = h.run();
        (status, h.get_model_status())
    }

    #[test]
    fn h16_eq_row_adds_on_sparsify_primal_dual_postsolve() {
        let (status, model_status) = sparsify_primal_dual_postsolve(true);
        assert_eq!(status, HighsStatus::kOk);
        assert_eq!(model_status, HighsModelStatus::kOptimal);
    }
}

#[cfg(test)]
mod x1b_psr_a_tests {
    use super::*;
    use crate::lp_data::highs_options::HighsOptions;

    #[test]
    fn x1b_psr_a_append_col_and_undo_resize() {
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(2, 2);
        stack.append_col_to_model();
        assert_eq!(stack.get_next_col_index(), 3);
        assert!(!stack.is_orig_col(2));
        assert_eq!(stack.get_orig_cols(), vec![0, 1]);

        let mut solution = HighsSolution::default();
        solution.value_valid = true;
        solution.dual_valid = true;
        solution.col_value = vec![1.0, 2.0, 3.0];
        solution.col_dual = vec![0.0; 3];
        solution.row_value = vec![0.0; 2];
        solution.row_dual = vec![0.0; 2];
        let mut basis = HighsBasis::default();
        let mut options = HighsOptions::default();
        options.mip_presolve_row_append = true;
        stack.undo(&options, &mut solution, &mut basis, -1, false);
        assert_eq!(solution.col_value, vec![1.0, 2.0]);
        assert_eq!(solution.col_dual.len(), 2);
        assert_eq!(stack.get_next_col_index(), 3);
    }

    #[test]
    fn x1b_psr_a_cut_row_types() {
        let mut stack = HighsPostsolveStack::default();
        stack.initialize_index_maps(2, 2);
        stack.append_cuts_to_model_x1b(2);
        stack.append_rows_to_model(1);
        assert_eq!(stack.get_orig_num_row(), 4);
        assert_eq!(stack.get_next_row_index(), 5);
        assert!(stack.is_orig_row(0) && stack.is_orig_row(1));
        assert!(stack.is_cut_row(2) && stack.is_cut_row(3));
        assert!(stack.is_appended_row(4));
        assert_eq!(stack.get_orig_rows(), vec![0, 1]);
        assert_eq!(stack.get_cut_rows(), vec![3, 2]);
        assert_eq!(stack.get_non_cut_rows(), vec![0, 1, 4]);

        stack.remove_cuts_from_model_x1b(2);
        assert_eq!(stack.get_orig_num_row(), 2);
        assert_eq!(stack.get_next_row_index(), 5);
        assert_eq!(stack.get_orig_rows_index(), &[0, 1, 4]);
        assert!(stack.is_appended_row(2));
        assert!(stack.get_cut_rows().is_empty());
    }
}
