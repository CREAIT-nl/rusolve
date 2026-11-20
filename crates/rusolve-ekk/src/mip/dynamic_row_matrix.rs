//! HiGHS source: mip/HighsDynamicRowMatrix.h, mip/HighsDynamicRowMatrix.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix (lines 17-102)
#[derive(Clone, Default, Debug)]
pub struct HighsDynamicRowMatrix {
    /// C++: `std::vector<std::pair<HighsInt, HighsInt>> ARrange_`
    pub a_rrange_: Vec<(HighsInt, HighsInt)>,
    /// C++: `std::vector<HighsInt> ARindex_`
    pub a_rindex_: Vec<HighsInt>,
    /// C++: `std::vector<double> ARvalue_`
    pub a_rvalue_: Vec<f64>,
    /// C++: `std::vector<HighsInt> ARrowindex_`
    pub a_rrowindex_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> AnextPos_`
    pub anext_pos_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> AprevPos_`
    pub aprev_pos_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> AnextNeg_`
    pub anext_neg_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> AprevNeg_`
    pub aprev_neg_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> AheadPos_`
    pub ahead_pos_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> AheadNeg_`
    pub ahead_neg_: Vec<HighsInt>,
    /// C++: `std::vector<uint8_t> colsLinked`
    pub cols_linked: Vec<u8>,
    /// C++: `std::set<std::pair<HighsInt, HighsInt>> freespaces_`
    pub freespaces_: BTreeSet<(HighsInt, HighsInt)>,
    /// C++: `std::vector<HighsInt> deletedrows_`
    pub deletedrows_: Vec<HighsInt>,
}

impl HighsDynamicRowMatrix {
    /// HiGHS: HighsDynamicRowMatrix.cpp::HighsDynamicRowMatrix::HighsDynamicRowMatrix (lines 15-18)
    /// C++: `HighsDynamicRowMatrix(HighsInt ncols)`
    pub fn new(ncols: HighsInt) -> Self {
        Self {
            ahead_pos_: vec![-1; ncols as usize],
            ahead_neg_: vec![-1; ncols as usize],
            ..Default::default()
        }
    }

    /// HiGHS: HighsDynamicRowMatrix.cpp::HighsDynamicRowMatrix::addRow (lines 20-98)
    /// C++: `HighsInt addRow(HighsInt * Rindex, double * Rvalue, HighsInt Rlen, bool linkCols = true)`
    /// C++ default arguments: `link_cols = true`
    pub fn add_row(
        &mut self,
        rindex: &mut [HighsInt],
        rvalue: &mut [f64],
        rlen: HighsInt,
        link_cols: bool,
    ) -> HighsInt {
        let rowindex: HighsInt;
        let start: HighsInt;
        let end: HighsInt;

        // insert the row in an existing empty space or append the values to the end
        // if no space that is large enough exists
        // C++: freespaces_.empty() || (it = freespaces_.lower_bound(make_pair(Rlen, -1))) == end()
        let it = if self.freespaces_.is_empty() {
            None
        } else {
            self.freespaces_.range((rlen, -1)..).next().copied()
        };
        match it {
            None => {
                start = self.a_rindex_.len() as HighsInt;
                end = start + rlen;

                let n = end as usize;
                self.a_rindex_.resize(n, 0);
                self.a_rvalue_.resize(n, 0.0);
                self.a_rrowindex_.resize(n, 0);
                self.aprev_pos_.resize(n, -1);
                self.anext_pos_.resize(n, -1);
                self.aprev_neg_.resize(n, -1);
                self.anext_neg_.resize(n, -1);
            }
            Some(freeslot) => {
                self.freespaces_.remove(&freeslot);

                start = freeslot.1;
                end = start + rlen;
                // if the space was not completely occupied, we register the remainder of
                // it again in the priority queue
                if freeslot.0 > rlen {
                    self.freespaces_.insert((freeslot.0 - rlen, end));
                }
            }
        }

        // register the range of values for this row with a reused or a new index
        if self.deletedrows_.is_empty() {
            rowindex = self.a_rrange_.len() as HighsInt;
            self.a_rrange_.push((start, end));
            self.cols_linked.push(link_cols as u8);
        } else {
            rowindex = *self.deletedrows_.last().expect("deletedrows_ is non-empty");
            self.deletedrows_.pop();
            self.a_rrange_[rowindex as usize].0 = start;
            self.a_rrange_[rowindex as usize].1 = end;
            self.cols_linked[rowindex as usize] = link_cols as u8;
        }

        // now add the nonzeros in the order sorted by the index value
        for i in start..end {
            let k = (i - start) as usize;
            self.a_rindex_[i as usize] = rindex[k];
            self.a_rvalue_[i as usize] = rvalue[k];
            self.a_rrowindex_[i as usize] = rowindex;
        }

        // link the row values to the columns
        if !link_cols {
            return rowindex;
        }

        for i in start..end {
            let col = self.a_rindex_[i as usize];

            if self.a_rvalue_[i as usize] > 0.0 {
                self.aprev_pos_[i as usize] = -1;
                let head = self.ahead_pos_[col as usize];
                self.ahead_pos_[col as usize] = i;
                self.anext_pos_[i as usize] = head;

                if head != -1 {
                    self.aprev_pos_[head as usize] = i;
                }
            } else {
                self.aprev_neg_[i as usize] = -1;
                let head = self.ahead_neg_[col as usize];
                self.ahead_neg_[col as usize] = i;
                self.anext_neg_[i as usize] = head;

                if head != -1 {
                    self.aprev_neg_[head as usize] = i;
                }
            }
        }

        rowindex
    }

    /// HiGHS: HighsDynamicRowMatrix.cpp::HighsDynamicRowMatrix::unlinkColumns (lines 100-143)
    /// C++: `void unlinkColumns(HighsInt rowindex)`
    pub fn unlink_columns(&mut self, rowindex: HighsInt) {
        if self.cols_linked[rowindex as usize] == 0 {
            return;
        }

        self.cols_linked[rowindex as usize] = 0;
        let start = self.a_rrange_[rowindex as usize].0;
        let end = self.a_rrange_[rowindex as usize].1;
        for i in start..end {
            let col = self.a_rindex_[i as usize];

            if self.a_rvalue_[i as usize] > 0.0 {
                let prev = self.aprev_pos_[i as usize];
                let next = self.anext_pos_[i as usize];

                if next != -1 {
                    debug_assert!(self.aprev_pos_[next as usize] == i);
                    self.aprev_pos_[next as usize] = prev;
                }

                if prev != -1 {
                    debug_assert!(self.anext_pos_[prev as usize] == i);
                    self.anext_pos_[prev as usize] = next;
                } else {
                    debug_assert!(self.ahead_pos_[col as usize] == i);
                    self.ahead_pos_[col as usize] = next;
                }
            } else {
                let prev = self.aprev_neg_[i as usize];
                let next = self.anext_neg_[i as usize];

                if next != -1 {
                    debug_assert!(self.aprev_neg_[next as usize] == i);
                    self.aprev_neg_[next as usize] = prev;
                }

                if prev != -1 {
                    debug_assert!(self.anext_neg_[prev as usize] == i);
                    self.anext_neg_[prev as usize] = next;
                } else {
                    debug_assert!(self.ahead_neg_[col as usize] == i);
                    self.ahead_neg_[col as usize] = next;
                }
            }
        }
    }

    /// HiGHS: HighsDynamicRowMatrix.cpp::HighsDynamicRowMatrix::removeRow (lines 147-199)
    /// C++: `void removeRow(HighsInt rowindex)`
    pub fn remove_row(&mut self, rowindex: HighsInt) {
        let start = self.a_rrange_[rowindex as usize].0;
        let end = self.a_rrange_[rowindex as usize].1;

        if self.cols_linked[rowindex as usize] != 0 {
            for i in start..end {
                let col = self.a_rindex_[i as usize];

                if self.a_rvalue_[i as usize] > 0.0 {
                    let prev = self.aprev_pos_[i as usize];
                    let next = self.anext_pos_[i as usize];

                    if next != -1 {
                        debug_assert!(self.aprev_pos_[next as usize] == i);
                        self.aprev_pos_[next as usize] = prev;
                    }

                    if prev != -1 {
                        debug_assert!(self.anext_pos_[prev as usize] == i);
                        self.anext_pos_[prev as usize] = next;
                    } else {
                        debug_assert!(self.ahead_pos_[col as usize] == i);
                        self.ahead_pos_[col as usize] = next;
                    }
                } else {
                    let prev = self.aprev_neg_[i as usize];
                    let next = self.anext_neg_[i as usize];

                    if next != -1 {
                        debug_assert!(self.aprev_neg_[next as usize] == i);
                        self.aprev_neg_[next as usize] = prev;
                    }

                    if prev != -1 {
                        debug_assert!(self.anext_neg_[prev as usize] == i);
                        self.anext_neg_[prev as usize] = next;
                    } else {
                        debug_assert!(self.ahead_neg_[col as usize] == i);
                        self.ahead_neg_[col as usize] = next;
                    }
                }
            }
        }

        // register the space of the deleted row and the index so that it can be
        // reused
        self.deletedrows_.push(rowindex);
        self.freespaces_.insert((end - start, start));

        // set the range to -1,-1 to indicate a deleted row
        self.a_rrange_[rowindex as usize].0 = -1;
        self.a_rrange_[rowindex as usize].1 = -1;
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::columnsLinked (lines 52-54)
    /// C++: `bool columnsLinked(HighsInt rowindex) const`
    pub fn columns_linked(&self, rowindex: HighsInt) -> bool {
        self.cols_linked[rowindex as usize] != 0
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::nonzeroCapacity (lines 66-66)
    /// C++: `std::size_t nonzeroCapacity() const`
    pub fn nonzero_capacity(&self) -> usize {
        self.a_rvalue_.len()
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::forEachPositiveColumnEntry (lines 71-79)
    /// C++: `void forEachPositiveColumnEntry() const`
    pub fn for_each_positive_column_entry(
        &self,
        col: HighsInt,
        f: &mut dyn FnMut(HighsInt, f64) -> bool,
    ) {
        let mut iter = self.ahead_pos_[col as usize];

        while iter != -1 {
            if !f(
                self.a_rrowindex_[iter as usize],
                self.a_rvalue_[iter as usize],
            ) {
                break;
            }
            iter = self.anext_pos_[iter as usize];
        }
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::forEachNegativeColumnEntry (lines 81-89)
    /// C++: `void forEachNegativeColumnEntry() const`
    pub fn for_each_negative_column_entry(
        &self,
        col: HighsInt,
        f: &mut dyn FnMut(HighsInt, f64) -> bool,
    ) {
        let mut iter = self.ahead_neg_[col as usize];

        while iter != -1 {
            if !f(
                self.a_rrowindex_[iter as usize],
                self.a_rvalue_[iter as usize],
            ) {
                break;
            }
            iter = self.anext_neg_[iter as usize];
        }
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::getNumRows (lines 91-91)
    /// C++: `HighsInt getNumRows() const`
    pub fn get_num_rows(&self) -> HighsInt {
        self.a_rrange_.len() as HighsInt
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::getNumDelRows (lines 93-93)
    /// C++: `HighsInt getNumDelRows() const`
    pub fn get_num_del_rows(&self) -> HighsInt {
        self.deletedrows_.len() as HighsInt
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::getRowStart (lines 95-95)
    /// C++: `HighsInt getRowStart(HighsInt row) const`
    pub fn get_row_start(&self, row: HighsInt) -> HighsInt {
        self.a_rrange_[row as usize].0
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::getRowEnd (lines 97-97)
    /// C++: `HighsInt getRowEnd(HighsInt row) const`
    pub fn get_row_end(&self, row: HighsInt) -> HighsInt {
        self.a_rrange_[row as usize].1
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::getARindex (lines 99-99)
    /// C++: `const HighsInt * getARindex() const`
    pub fn get_a_rindex(&self) -> &[HighsInt] {
        &self.a_rindex_
    }

    /// HiGHS: HighsDynamicRowMatrix.h::HighsDynamicRowMatrix::getARvalue (lines 101-101)
    /// C++: `const double * getARvalue() const`
    pub fn get_a_rvalue(&self) -> &[f64] {
        &self.a_rvalue_
    }
}
