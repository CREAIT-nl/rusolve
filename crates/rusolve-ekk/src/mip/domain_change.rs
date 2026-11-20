//! HiGHS source: mip/HighsDomainChange.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsDomainChange.h::HighsDomainChange (lines 16-39)
#[derive(Clone, Default, Debug)]
pub struct HighsDomainChange {
    /// C++: `double boundval`
    pub boundval: f64,
    /// C++: `HighsInt column`
    pub column: HighsInt,
    /// C++: `HighsBoundType boundtype`
    pub boundtype: HighsBoundType,
}

/// `std::set`/`std::sort` ordering = C++ `operator<` (equivalence = neither is less).
impl PartialEq for HighsDomainChange {
    fn eq(&self, other: &Self) -> bool {
        self.op_eq(other)
    }
}
impl Eq for HighsDomainChange {}
impl PartialOrd for HighsDomainChange {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HighsDomainChange {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.op_lt(other) {
            std::cmp::Ordering::Less
        } else if other.op_lt(self) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    }
}

/// HiGHS: HighsDomainChange.h::HighsSubstitution (lines 41-46)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct HighsSubstitution {
    /// C++: `HighsInt substcol`
    pub substcol: HighsInt,
    /// C++: `HighsInt staycol`
    pub staycol: HighsInt,
    /// C++: `double scale`
    pub scale: f64,
    /// C++: `double offset`
    pub offset: f64,
}

/// HiGHS: HighsDomainChange.h::HighsBoundType (enum class)
#[repr(i32)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub enum HighsBoundType {
    #[default]
    kLower = 0,
    kUpper = 1,
}

impl HighsDomainChange {
    /// HiGHS: HighsDomainChange.h::HighsDomainChange::operator< (lines 21-28)
    /// C++: `bool operator<(const HighsDomainChange & other) const`
    pub fn op_lt(&self, other: &HighsDomainChange) -> bool {
        if self.column < other.column {
            return true;
        }
        if other.column < self.column {
            return false;
        }
        if (self.boundtype as HighsInt) < (other.boundtype as HighsInt) {
            return true;
        }
        if (other.boundtype as HighsInt) < (self.boundtype as HighsInt) {
            return false;
        }
        if self.boundval < other.boundval {
            return true;
        }
        false
    }

    /// HiGHS: HighsDomainChange.h::HighsDomainChange::operator== (lines 30-33)
    /// C++: `bool operator==(const HighsDomainChange & other) const`
    pub fn op_eq(&self, other: &HighsDomainChange) -> bool {
        self.boundtype == other.boundtype
            && self.column == other.column
            && self.boundval == other.boundval
    }

    /// HiGHS: HighsDomainChange.h::HighsDomainChange::operator!= (lines 35-38)
    /// C++: `bool operator!=(const HighsDomainChange & other) const`
    pub fn op_ne(&self, other: &HighsDomainChange) -> bool {
        self.boundtype != other.boundtype
            || self.column != other.column
            || self.boundval != other.boundval
    }
}
