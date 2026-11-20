//! HiGHS source: util/HighsDisjointSets.h (lines 1-107)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsDisjointSets.h::HighsDisjointSets<kMinimalRepresentative> (lines 27-105)
///
/// `template <bool kMinimalRepresentative = false>` -> const generic with the
/// same default (`HighsDisjointSets` = `HighsDisjointSets<false>`).
#[derive(Clone, Debug, Default)]
pub struct HighsDisjointSets<const K_MINIMAL_REPRESENTATIVE: bool = false> {
    /// C++: `std::vector<HighsInt> sizes`
    pub sizes: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> sets`
    pub sets: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> linkCompressionStack`
    pub link_compression_stack: Vec<HighsInt>,
}

impl<const K_MINIMAL_REPRESENTATIVE: bool> HighsDisjointSets<K_MINIMAL_REPRESENTATIVE> {
    /// HiGHS: HighsDisjointSets.h::HighsDisjointSets::HighsDisjointSets (line 36)
    /// C++: `HighsDisjointSets(HighsInt numItems)`
    pub fn new(num_items: HighsInt) -> Self {
        let mut s = Self::default();
        s.reset(num_items);
        s
    }

    /// HiGHS: HighsDisjointSets.h::HighsDisjointSets::reset (lines 38-42)
    pub fn reset(&mut self, num_items: HighsInt) {
        let n = num_items as usize;
        // sizes.assign(numItems, 1)
        self.sizes.clear();
        self.sizes.resize(n, 1);
        // sets.resize(numItems); std::iota(sets.begin(), sets.end(), 0);
        self.sets.resize(n, 0);
        for (i, s) in self.sets.iter_mut().enumerate() {
            *s = i as HighsInt;
        }
    }

    /// HiGHS: HighsDisjointSets.h::HighsDisjointSets::getSet (lines 44-66)
    pub fn get_set(&mut self, item: HighsInt) -> HighsInt {
        let mut item = item;
        debug_assert!(item >= 0 && (item as usize) < self.sets.len());
        let mut repr = self.sets[item as usize];
        debug_assert!(repr >= 0 && (repr as usize) < self.sets.len());

        if repr != self.sets[repr as usize] {
            // do { push(item); item = repr; repr = sets[repr]; } while (repr != sets[repr]);
            loop {
                self.link_compression_stack.push(item);
                item = repr;
                repr = self.sets[repr as usize];
                if repr == self.sets[repr as usize] {
                    break;
                }
            }

            // do { i = back(); pop_back(); sets[i] = repr; } while (!empty());
            loop {
                let i = self.link_compression_stack.pop().unwrap();
                self.sets[i as usize] = repr;
                if self.link_compression_stack.is_empty() {
                    break;
                }
            }

            self.sets[item as usize] = repr;
        }

        repr
    }

    /// HiGHS: HighsDisjointSets.h::HighsDisjointSets::getSetSize (lines 68-71)
    pub fn get_set_size(&self, set: HighsInt) -> HighsInt {
        debug_assert!(self.sets[set as usize] == set);
        self.sizes[set as usize]
    }

    /// HiGHS: HighsDisjointSets.h::HighsDisjointSets::merge (lines 73-104)
    pub fn merge(&mut self, item1: HighsInt, item2: HighsInt) {
        debug_assert!(item1 >= 0 && (item1 as usize) < self.sets.len());
        debug_assert!(item2 >= 0 && (item2 as usize) < self.sets.len());

        let repr1 = self.get_set(item1);
        debug_assert!(self.sets[repr1 as usize] == repr1);
        debug_assert!(repr1 >= 0 && (repr1 as usize) < self.sets.len());

        let repr2 = self.get_set(item2);
        debug_assert!(self.sets[repr2 as usize] == repr2);
        debug_assert!(repr2 >= 0 && (repr2 as usize) < self.sets.len());
        debug_assert!(self.sizes.len() == self.sets.len());

        if repr1 == repr2 {
            return;
        }

        let r1 = repr1 as usize;
        let r2 = repr2 as usize;
        if K_MINIMAL_REPRESENTATIVE {
            if repr2 > repr1 {
                self.sets[r2] = repr1;
                let s2 = self.sizes[r2];
                self.sizes[r1] += s2;
            } else {
                self.sets[r1] = repr2;
                let s1 = self.sizes[r1];
                self.sizes[r2] += s1;
            }
        } else if self.sizes[r1] > self.sizes[r2] {
            self.sets[r2] = repr1;
            let s2 = self.sizes[r2];
            self.sizes[r1] += s2;
        } else {
            self.sets[r1] = repr2;
            let s1 = self.sizes[r1];
            self.sizes[r2] += s1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_and_find() {
        let mut d = HighsDisjointSets::<false>::new(5);
        d.merge(0, 1);
        d.merge(2, 3);
        d.merge(1, 3);
        let r = d.get_set(0);
        assert_eq!(r, d.get_set(2));
        assert_eq!(d.get_set_size(r), 4);
        let single = d.get_set(4);
        assert_eq!(d.get_set_size(single), 1);
    }

    #[test]
    fn minimal_representative_is_smallest_index() {
        let mut d = HighsDisjointSets::<true>::new(4);
        d.merge(3, 1);
        assert_eq!(d.get_set(3), 1);
        d.merge(0, 3);
        assert_eq!(d.get_set(0), 0);
        assert_eq!(d.get_set(1), 0);
    }
}
