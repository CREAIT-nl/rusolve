//! HiGHS source: util/HighsMatrixSlice.h (HighsSliceNonzero, HighsEmptySlice,
//! HighsTripletListSlice, HighsTripletTreeSlicePreOrder,
//! HighsTripletTreeSliceInOrder, HighsTripletPositionSlice)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Complete in the skeleton (no fill-in needed).
//!
//! # Design
//!
//! C++ slices hold raw pointers into HPresolve's triplet storage (`Avalue`,
//! `Arow`, `Acol`, `Anext`, `ARleft`, `ARright`, `rowpositions`) and their
//! iterators read those arrays *lazily*: `operator++` follows the link of the
//! current node only when the range-for advances, i.e. after the loop body
//! ran. Rust cannot hold those borrows across a body that mutates HPresolve,
//! so a slice here is a small **cursor** (indices only, no borrow) whose
//! [`HighsTripletSlice::next_nz`] takes the storage as an argument and does
//! exactly what C++ `++it; *it` does at that moment:
//!
//! * the first call performs `begin()` (for the in-order tree slice this is
//!   the descent to the leftmost node), later calls perform `operator++`,
//!   reading the links *now* (after the previous body);
//! * it returns the nonzero (`index`, `value`, `position`) of the new current
//!   node, or `None` at `end()`.
//!
//! The storage is anything implementing [`HighsTripletStorage`]; HPresolve
//! implements it. Two ways to loop:
//!
//! ```ignore
//! // Body only reads `self`: borrowing iterator (C++ range-for).
//! for nz in self.get_column_vector(col).iter(self) { ... nz.index() ... }
//!
//! // Body mutates `self`: keep the cursor, pass `self` at each step.
//! let mut it = self.get_row_vector(row);
//! while let Some(nz) = it.next_nz(self) { ... self.unlink(nz.position()) ... }
//! ```
//!
//! Caveat (the one semantic difference): `HighsSliceNonzero` holds *copies*
//! of `index`/`value` taken at `next_nz`, while C++ `nonzero.value()`
//! dereferences `Avalue[pos]` every time. If a loop body writes
//! `Avalue[nz.position()]` and then reads `nonzero.value()` again, read
//! `self.a_value[nz.position() as usize]` instead.
//!
//! `iterator.position()` (C++, used by `getRowPositions`) is
//! [`HighsSliceNonzero::position`]. The compressed/indexed slices are not used
//! by LP presolve and are omitted. `HighsEmptySlice()` is [`HighsEmptySlice`].
//!
//! Postsolve-stack methods that take a `const HighsMatrixSlice<F>&` take
//! `impl IntoIterator<Item = HighsSliceNonzero>`: pass `slice.iter(self)` or
//! `HighsEmptySlice`.

use crate::prelude::*;

/// HiGHS: HighsMatrixSlice.h::HighsSliceNonzero
///
/// `position` is the C++ iterator's `position()` (the node / nonzero slot).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HighsSliceNonzero {
    index_: HighsInt,
    value_: f64,
    position_: HighsInt,
}

impl HighsSliceNonzero {
    pub fn new(index: HighsInt, value: f64, position: HighsInt) -> Self {
        Self {
            index_: index,
            value_: value,
            position_: position,
        }
    }
    /// HiGHS: HighsSliceNonzero::index
    pub fn index(&self) -> HighsInt {
        self.index_
    }
    /// HiGHS: HighsSliceNonzero::value
    pub fn value(&self) -> f64 {
        self.value_
    }
    /// HiGHS: `iterator::position()` of the iterator that produced this nonzero.
    pub fn position(&self) -> HighsInt {
        self.position_
    }
}

/// Borrowed view of HPresolve's triplet storage (C++ raw pointers).
#[derive(Clone, Copy)]
pub struct HighsTripletArrays<'a> {
    pub a_value: &'a [f64],
    pub a_row: &'a [HighsInt],
    pub a_col: &'a [HighsInt],
    pub a_next: &'a [HighsInt],
    pub a_r_left: &'a [HighsInt],
    pub a_r_right: &'a [HighsInt],
    pub row_positions: &'a [HighsInt],
}

/// Owner of triplet storage that slices iterate over (HPresolve).
pub trait HighsTripletStorage {
    fn triplet_arrays(&self) -> HighsTripletArrays<'_>;
}

/// A lazy cursor over one row or column of triplet storage.
pub trait HighsTripletSlice: Sized {
    /// C++ `begin()` on the first call, `operator++` afterwards, then `*it`;
    /// `None` at `end()`.
    fn next_nz<S: HighsTripletStorage + ?Sized>(&mut self, s: &S) -> Option<HighsSliceNonzero>;

    /// Borrowing iterator for loop bodies that do not mutate the storage.
    fn iter<S: HighsTripletStorage + ?Sized>(self, s: &S) -> HighsSliceIter<'_, Self, S> {
        HighsSliceIter {
            slice: self,
            storage: s,
        }
    }
}

/// Iterator adaptor returned by [`HighsTripletSlice::iter`].
pub struct HighsSliceIter<'a, T, S: ?Sized> {
    slice: T,
    storage: &'a S,
}

impl<T: HighsTripletSlice, S: HighsTripletStorage + ?Sized> Iterator for HighsSliceIter<'_, T, S> {
    type Item = HighsSliceNonzero;
    fn next(&mut self) -> Option<HighsSliceNonzero> {
        self.slice.next_nz(self.storage)
    }
}

/// HiGHS: HighsMatrixSlice.h::HighsEmptySlice (`begin() == end()`).
#[derive(Clone, Copy, Debug, Default)]
pub struct HighsEmptySlice;

impl IntoIterator for HighsEmptySlice {
    type Item = HighsSliceNonzero;
    type IntoIter = std::iter::Empty<HighsSliceNonzero>;
    fn into_iter(self) -> Self::IntoIter {
        std::iter::empty()
    }
}

/// HiGHS: HighsMatrixSlice.h::HighsMatrixSlice<HighsTripletListSlice>
/// (column linked list: index = `Arow`, value = `Avalue`, next = `Anext`).
#[derive(Clone, Debug)]
pub struct HighsTripletListSlice {
    head: HighsInt,
    current: HighsInt,
    started: bool,
}

impl HighsTripletListSlice {
    pub fn new(head: HighsInt) -> Self {
        Self {
            head,
            current: -1,
            started: false,
        }
    }
}

impl HighsTripletSlice for HighsTripletListSlice {
    fn next_nz<S: HighsTripletStorage + ?Sized>(&mut self, s: &S) -> Option<HighsSliceNonzero> {
        let a = s.triplet_arrays();
        if !self.started {
            self.started = true;
            self.current = self.head;
        } else {
            if self.current == -1 {
                return None;
            }
            self.current = a.a_next[self.current as usize];
        }
        if self.current == -1 {
            return None;
        }
        let pos = self.current as usize;
        Some(HighsSliceNonzero::new(
            a.a_row[pos],
            a.a_value[pos],
            self.current,
        ))
    }
}

/// Traversal stack of the tree slices: C++ `std::vector<HighsInt> stack` with
/// `stack.reserve(16)` at `begin()`.
///
/// PORT NOTE (perf): the first [`SLICE_STACK_INLINE`] entries live inline, so a
/// row walk allocates nothing unless the splay tree is deeper than that (the
/// C++ iterator mallocs its reserve(16) on every `begin()`). Push/pop order is
/// the vector's, so the visit order is unchanged.
pub const SLICE_STACK_INLINE: usize = 32;

#[derive(Clone, Debug)]
pub struct SliceStack {
    small: [HighsInt; SLICE_STACK_INLINE],
    len: usize,
    big: Vec<HighsInt>,
}

impl Default for SliceStack {
    fn default() -> Self {
        Self::new()
    }
}

impl SliceStack {
    #[inline]
    pub fn new() -> Self {
        Self {
            small: [0; SLICE_STACK_INLINE],
            len: 0,
            big: Vec::new(),
        }
    }

    #[inline(always)]
    pub fn push(&mut self, v: HighsInt) {
        if self.len < SLICE_STACK_INLINE {
            self.small[self.len] = v;
        } else {
            self.big.push(v);
        }
        self.len += 1;
    }

    /// `back(); pop_back()`; `None` when empty.
    #[inline(always)]
    pub fn pop(&mut self) -> Option<HighsInt> {
        if self.len == 0 {
            return None;
        }
        self.len -= 1;
        if self.len < SLICE_STACK_INLINE {
            Some(self.small[self.len])
        } else {
            self.big.pop()
        }
    }
}

/// HiGHS: HighsMatrixSlice.h::HighsMatrixSlice<HighsTripletTreeSlicePreOrder>
/// (row splay tree, pre-order: index = `Acol`, value = `Avalue`).
#[derive(Clone, Debug)]
pub struct HighsTripletTreeSlicePreOrder {
    root: HighsInt,
    current: HighsInt,
    stack: SliceStack,
    started: bool,
}

impl HighsTripletTreeSlicePreOrder {
    pub fn new(root: HighsInt) -> Self {
        Self {
            root,
            current: -1,
            stack: SliceStack::new(),
            started: false,
        }
    }
}

impl HighsTripletSlice for HighsTripletTreeSlicePreOrder {
    fn next_nz<S: HighsTripletStorage + ?Sized>(&mut self, s: &S) -> Option<HighsSliceNonzero> {
        let a = s.triplet_arrays();
        if !self.started {
            self.started = true;
            if self.root == -1 {
                self.current = -1;
                return None;
            }
            self.current = self.root;
            self.stack.push(-1);
        } else {
            if self.current == -1 {
                return None;
            }
            let cur = self.current as usize;
            if a.a_r_left[cur] != -1 {
                if a.a_r_right[cur] != -1 {
                    self.stack.push(a.a_r_right[cur]);
                }
                self.current = a.a_r_left[cur];
            } else if a.a_r_right[cur] != -1 {
                self.current = a.a_r_right[cur];
            } else {
                self.current = self.stack.pop().unwrap_or(-1);
            }
        }
        if self.current == -1 {
            return None;
        }
        let pos = self.current as usize;
        Some(HighsSliceNonzero::new(
            a.a_col[pos],
            a.a_value[pos],
            self.current,
        ))
    }
}

/// HiGHS: HighsMatrixSlice.h::HighsMatrixSlice<HighsTripletTreeSliceInOrder>
/// (row splay tree, in-order = increasing column: index = `Acol`, value = `Avalue`).
#[derive(Clone, Debug)]
pub struct HighsTripletTreeSliceInOrder {
    root: HighsInt,
    current: HighsInt,
    stack: SliceStack,
    started: bool,
}

impl HighsTripletTreeSliceInOrder {
    pub fn new(root: HighsInt) -> Self {
        Self {
            root,
            current: -1,
            stack: SliceStack::new(),
            started: false,
        }
    }
}

impl HighsTripletSlice for HighsTripletTreeSliceInOrder {
    fn next_nz<S: HighsTripletStorage + ?Sized>(&mut self, s: &S) -> Option<HighsSliceNonzero> {
        let a = s.triplet_arrays();
        if !self.started {
            self.started = true;
            if self.root == -1 {
                self.current = -1;
                return None;
            }
            self.current = self.root;
            self.stack.push(-1);
            while a.a_r_left[self.current as usize] != -1 {
                self.stack.push(self.current);
                self.current = a.a_r_left[self.current as usize];
            }
        } else {
            if self.current == -1 {
                return None;
            }
            let cur = self.current as usize;
            if a.a_r_right[cur] != -1 {
                self.current = a.a_r_right[cur];
                while a.a_r_left[self.current as usize] != -1 {
                    self.stack.push(self.current);
                    self.current = a.a_r_left[self.current as usize];
                }
            } else {
                self.current = self.stack.pop().unwrap_or(-1);
            }
        }
        if self.current == -1 {
            return None;
        }
        let pos = self.current as usize;
        Some(HighsSliceNonzero::new(
            a.a_col[pos],
            a.a_value[pos],
            self.current,
        ))
    }
}

/// HiGHS: HighsMatrixSlice.h::HighsMatrixSlice<HighsTripletPositionSlice>
/// (the stored row: positions from `rowpositions`, index = `Acol`,
/// value = `Avalue`). `len` is captured at creation, as in C++.
#[derive(Clone, Debug)]
pub struct HighsTripletPositionSlice {
    len: usize,
    i: usize,
    started: bool,
}

impl HighsTripletPositionSlice {
    pub fn new(len: usize) -> Self {
        Self {
            len,
            i: 0,
            started: false,
        }
    }
}

impl HighsTripletSlice for HighsTripletPositionSlice {
    fn next_nz<S: HighsTripletStorage + ?Sized>(&mut self, s: &S) -> Option<HighsSliceNonzero> {
        let a = s.triplet_arrays();
        if !self.started {
            self.started = true;
        } else if self.i < self.len {
            self.i += 1;
        }
        if self.i >= self.len {
            return None;
        }
        let node = a.row_positions[self.i];
        let pos = node as usize;
        Some(HighsSliceNonzero::new(a.a_col[pos], a.a_value[pos], node))
    }
}
