//! HiGHS source: qpsolver/snippets.hpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Generic helpers; reached as `crate::qpsolver::snippets::contains`
//! etc. (the names are generic, so they are NOT re-exported by `qpsolver::*`).

use crate::prelude::*;

/// HiGHS: snippets.hpp::contains (lines 14-17)
pub fn contains<T: PartialEq>(vec: &[T], element: &T) -> bool {
    vec.contains(element)
}

/// HiGHS: snippets.hpp::remove (lines 19-24)
///
/// Returns what C++ returns: `rem2 != vec.end()` after the erase, i.e.
/// always `false` (the erase returns the new end).
pub fn remove<T: PartialEq>(vec: &mut Vec<T>, element: &T) -> bool {
    vec.retain(|x| *x != *element);
    false
}

/// HiGHS: snippets.hpp::indexof (lines 26-34)
pub fn indexof<T: PartialEq>(vec: &[T], element: &T) -> HighsInt {
    for (i, x) in vec.iter().enumerate() {
        if *x == *element {
            return i as HighsInt;
        }
    }
    -1
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use crate::qpsolver::*;

    #[test]
    fn find_and_remove() {
        let mut v: Vec<HighsInt> = vec![4, 7, 9];
        assert!(contains(&v, &7));
        assert_eq!(indexof(&v, &9), 2);
        assert_eq!(indexof(&v, &5), -1);
        // C++ remove (snippets.hpp:19-24) returns rem2 != vec.end() which is
        // always false after erase; the skeleton test asserted true. C++ wins.
        assert!(!remove(&mut v, &7));
        assert_eq!(v, vec![4, 9]);
        assert!(!remove(&mut v, &7));
    }
}
