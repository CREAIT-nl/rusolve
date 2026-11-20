//! HiGHS source: util/HighsHashTree.h (lines 1-1461)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//!
//! `HighsHashTree<K, V>` is a hash array mapped trie (HAMT). Its iteration
//! order (`for_each`) and `find_common` result depend on the node layout,
//! and HiGHS iterates it in `HighsImplications` (vub/vlb maps) and
//! `HighsCliqueTable` (inverted clique lists), so the port must reproduce the
//! **same tree shape**: same hash chunking (`kBitsPerLevel = 6`,
//! `kMaxDepth = 9`), same leaf size classes (`kMinLeafSize = 6`,
//! `kLeafBurstThreshold = 54`, `InnerLeaf<1..4>` capacities), same
//! occupation-bitmap child order, same list-leaf collision chains, same
//! insertion position inside leaves (sorted by the 16-bit hash chunk).
//!
//! C++ stores `NodePtr` as a tagged `uintptr_t` (type in the low 3 bits)
//! and allocates branch nodes with trailing children (`children()` =
//! `this + 1`). In Rust the tag + pointer is the enum [`HighsHashTreeNode`],
//! the trailing children a `Vec`, and the list-leaf `ListNode* next` chain
//! an `Option<Box<..>>`. No `unsafe`.
//!
//! `HighsHashTree<K>` (C++ `V = HighsHashTableEntry`'s `void` default) is
//! `HighsHashTree<K, ()>`.

use crate::prelude::*;

/// HiGHS: HighsHashTree.h::HighsHashTree::Type (lines 24-32)
pub const kEmpty: i32 = 0;
pub const kListLeaf: i32 = 1;
pub const kInnerLeafSizeClass1: i32 = 2;
pub const kInnerLeafSizeClass2: i32 = 3;
pub const kInnerLeafSizeClass3: i32 = 4;
pub const kInnerLeafSizeClass4: i32 = 5;
pub const kBranchNode: i32 = 6;
/// HiGHS: HighsHashTree.h::HighsHashTree::Constants (lines 34-44)
pub const kBitsPerLevel: i32 = 6;
pub const kBranchFactor: i32 = 1 << kBitsPerLevel;
pub const kMaxDepth: i32 = 9;
pub const kMinLeafSize: i32 = 6;
pub const kLeafBurstThreshold: i32 = 54;

/// HiGHS: HighsHashTree.h::HighsHashTree::Occupation (lines 69-91)
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HighsHashTree_Occupation {
    /// C++: `uint64_t occupation`
    pub occupation: u64,
}

impl HighsHashTree_Occupation {
    /// HiGHS: HighsHashTree.h::Occupation::set (line 76)
    pub fn set(&mut self, pos: u8) {
        self.occupation |= 1u64 << pos;
    }
    /// HiGHS: HighsHashTree.h::Occupation::flip (line 78)
    pub fn flip(&mut self, pos: u8) {
        self.occupation ^= 1u64 << pos;
    }
    /// HiGHS: HighsHashTree.h::Occupation::test (line 80)
    pub fn test(&self, pos: u8) -> bool {
        (self.occupation & (1u64 << pos)) != 0
    }
    /// HiGHS: HighsHashTree.h::Occupation::num_set_until (lines 82-84)
    pub fn num_set_until(&self, pos: u8) -> i32 {
        (self.occupation >> pos).count_ones() as i32
    }
    /// HiGHS: HighsHashTree.h::Occupation::num_set_after (lines 86-88)
    pub fn num_set_after(&self, pos: u8) -> i32 {
        (self.occupation << (63 - pos as u32)).count_ones() as i32
    }
    /// HiGHS: HighsHashTree.h::Occupation::num_set (line 90)
    pub fn num_set(&self) -> i32 {
        self.occupation.count_ones() as i32
    }
}

/// HiGHS: HighsHashTree.h::HighsHashTree::InnerLeaf<kSizeClass> (lines 99-274)
///
/// `kSizeClass` (1..=4) is carried by the capacity `CAP` = `capacity()` =
/// `kMinLeafSize + (kSizeClass - 1) * (kLeafBurstThreshold - kMinLeafSize) / 3`, with
/// `HCAP = CAP + 1` (stable Rust has no `CAP + 1` in an array length). As in C++,
/// `std::array<uint64_t, capacity() + 1> hashes` and `std::array<Entry, capacity()>
/// entries` are inline arrays, so a leaf is a single allocation (its `Box`). The four
/// classes are [`HighsHashTree_InnerLeaf1`] .. [`HighsHashTree_InnerLeaf4`].
#[derive(Clone, Debug)]
pub struct HighsHashTree_InnerLeaf<K, V, const CAP: usize, const HCAP: usize> {
    /// C++: `Occupation occupation`
    pub occupation: HighsHashTree_Occupation,
    /// C++: `int size`
    pub size: i32,
    /// C++: `std::array<uint64_t, capacity() + 1> hashes`
    pub hashes: [u64; HCAP],
    /// C++: `std::array<Entry, capacity()> entries`
    pub entries: [HighsHashTableEntry<K, V>; CAP],
}

/// C++ `InnerLeaf<1>` (capacity 6)
pub type HighsHashTree_InnerLeaf1<K, V> = HighsHashTree_InnerLeaf<K, V, 6, 7>;
/// C++ `InnerLeaf<2>` (capacity 22)
pub type HighsHashTree_InnerLeaf2<K, V> = HighsHashTree_InnerLeaf<K, V, 22, 23>;
/// C++ `InnerLeaf<3>` (capacity 38)
pub type HighsHashTree_InnerLeaf3<K, V> = HighsHashTree_InnerLeaf<K, V, 38, 39>;
/// C++ `InnerLeaf<4>` (capacity 54)
pub type HighsHashTree_InnerLeaf4<K, V> = HighsHashTree_InnerLeaf<K, V, 54, 55>;

/// HiGHS: HighsHashTree.h::HighsHashTree::InnerLeaf::capacity (lines 101-103) for a
/// runtime `kSizeClass`.
pub const fn inner_leaf_capacity(size_class: i32) -> i32 {
    kMinLeafSize + (size_class - 1) * (kLeafBurstThreshold - kMinLeafSize) / 3
}

/// HiGHS: HighsHashTree.h::HighsHashTree::ListNode (lines 276-281)
#[derive(Clone, Debug, Default)]
pub struct HighsHashTree_ListNode<K, V> {
    /// C++: `ListNode* next`
    pub next: Option<Box<HighsHashTree_ListNode<K, V>>>,
    /// C++: `HighsHashTableEntry<K, V> entry`
    pub entry: HighsHashTableEntry<K, V>,
}

/// HiGHS: HighsHashTree.h::HighsHashTree::ListLeaf (lines 282-288)
#[derive(Clone, Debug, Default)]
pub struct HighsHashTree_ListLeaf<K, V> {
    /// C++: `ListNode first`
    pub first: HighsHashTree_ListNode<K, V>,
    /// C++: `int count`
    pub count: i32,
}

/// HiGHS: HighsHashTree.h::HighsHashTree::BranchNode (lines 392-413)
/// (children stored after the node in C++; a `Vec` here, `occupation.num_set()` long).
#[derive(Clone, Debug, Default)]
pub struct HighsHashTree_BranchNode<K, V> {
    /// C++: `Occupation occupation`
    pub occupation: HighsHashTree_Occupation,
    /// C++: `NodePtr children()[numChilds]`
    pub children: Vec<HighsHashTreeNode<K, V>>,
}

/// HiGHS: HighsHashTree.h::HighsHashTree::NodePtr (lines 292-390): tagged pointer as an enum.
#[derive(Clone, Debug, Default)]
pub enum HighsHashTreeNode<K, V> {
    /// `kEmpty`
    #[default]
    Empty,
    /// `kListLeaf`
    ListLeaf(Box<HighsHashTree_ListLeaf<K, V>>),
    /// `kInnerLeafSizeClass1`
    InnerLeaf1(Box<HighsHashTree_InnerLeaf1<K, V>>),
    /// `kInnerLeafSizeClass2`
    InnerLeaf2(Box<HighsHashTree_InnerLeaf2<K, V>>),
    /// `kInnerLeafSizeClass3`
    InnerLeaf3(Box<HighsHashTree_InnerLeaf3<K, V>>),
    /// `kInnerLeafSizeClass4`
    InnerLeaf4(Box<HighsHashTree_InnerLeaf4<K, V>>),
    /// `kBranchNode`
    Branch(Box<HighsHashTree_BranchNode<K, V>>),
}

/// Evaluates `$body` with `$leaf` bound to the boxed leaf of an inner-leaf node (a
/// different type per size class, so one arm per class, like the C++ `switch` with a
/// `case` per `kInnerLeafSizeClass*` calling the same template). `$other` is the
/// expression for every other node type.
macro_rules! with_inner_leaf {
    ($node:expr, $leaf:ident => $body:expr, _ => $other:expr) => {
        match $node {
            HighsHashTreeNode::InnerLeaf1($leaf) => $body,
            HighsHashTreeNode::InnerLeaf2($leaf) => $body,
            HighsHashTreeNode::InnerLeaf3($leaf) => $body,
            HighsHashTreeNode::InnerLeaf4($leaf) => $body,
            _ => $other,
        }
    };
}

impl<K, V> HighsHashTreeNode<K, V> {
    /// HiGHS: HighsHashTree.h::NodePtr::getType (the `k*` type tag)
    pub fn get_type(&self) -> i32 {
        match self {
            HighsHashTreeNode::Empty => kEmpty,
            HighsHashTreeNode::ListLeaf(_) => kListLeaf,
            HighsHashTreeNode::InnerLeaf1(_) => kInnerLeafSizeClass1,
            HighsHashTreeNode::InnerLeaf2(_) => kInnerLeafSizeClass2,
            HighsHashTreeNode::InnerLeaf3(_) => kInnerLeafSizeClass3,
            HighsHashTreeNode::InnerLeaf4(_) => kInnerLeafSizeClass4,
            HighsHashTreeNode::Branch(_) => kBranchNode,
        }
    }

    /// HiGHS: HighsHashTree.h::NodePtr::numEntriesEstimate (lines 327-345)
    fn num_entries_estimate(&self) -> i32 {
        match self {
            HighsHashTreeNode::Empty => 0,
            HighsHashTreeNode::ListLeaf(_) => 1,
            HighsHashTreeNode::InnerLeaf1(_) => inner_leaf_capacity(1),
            HighsHashTreeNode::InnerLeaf2(_) => inner_leaf_capacity(2),
            HighsHashTreeNode::InnerLeaf3(_) => inner_leaf_capacity(3),
            HighsHashTreeNode::InnerLeaf4(_) => inner_leaf_capacity(4),
            HighsHashTreeNode::Branch(_) => kBranchFactor,
        }
    }

    /// HiGHS: HighsHashTree.h::NodePtr::numEntries (lines 347-365)
    fn num_entries(&self) -> i32 {
        match self {
            HighsHashTreeNode::Empty => 0,
            HighsHashTreeNode::ListLeaf(l) => l.count,
            HighsHashTreeNode::InnerLeaf1(l) => l.size,
            HighsHashTreeNode::InnerLeaf2(l) => l.size,
            HighsHashTreeNode::InnerLeaf3(l) => l.size,
            HighsHashTreeNode::InnerLeaf4(l) => l.size,
            HighsHashTreeNode::Branch(_) => kBranchFactor,
        }
    }
}

/// HiGHS: HighsHashTree.h::HighsHashTree<K, V> (lines 18-1459)
#[derive(Clone, Debug, Default)]
pub struct HighsHashTree<K, V = ()> {
    /// C++: `NodePtr root`
    pub root: HighsHashTreeNode<K, V>,
}

/// Private mirrors of the C++ `InnerLeaf<kSizeClass>` member functions (lines 99-274),
/// instantiated per size class as in C++.
impl<K: HighsHashKey + Copy, V: Copy + Default, const CAP: usize, const HCAP: usize>
    HighsHashTree_InnerLeaf<K, V, CAP, HCAP>
{
    /// HiGHS: HighsHashTree.h::InnerLeaf::capacity (lines 101-103)
    #[inline]
    const fn capacity() -> i32 {
        CAP as i32
    }

    /// HiGHS: HighsHashTree.h::InnerLeaf::InnerLeaf() (lines 116-117): occupation 0, size 0,
    /// hashes[0] = 0 (C++ leaves the other slots uninitialised; zero / default here).
    fn new_boxed() -> Box<Self> {
        debug_assert_eq!(HCAP, CAP + 1);
        Box::new(HighsHashTree_InnerLeaf {
            occupation: HighsHashTree_Occupation::default(),
            size: 0,
            hashes: [0; HCAP],
            entries: [HighsHashTableEntry::default(); CAP],
        })
    }

    /// HiGHS: HighsHashTree.h::InnerLeaf::InnerLeaf(InnerLeaf<kOtherSize>&&) (lines 118-127),
    /// the C++ template move constructor from another size class.
    fn from_moved_leaf<const OCAP: usize, const OHCAP: usize>(
        other: &HighsHashTree_InnerLeaf<K, V, OCAP, OHCAP>,
    ) -> Box<Self> {
        let mut leaf = Self::new_boxed();
        debug_assert!(other.size <= Self::capacity());
        leaf.occupation = other.occupation;
        leaf.size = other.size;
        let n = other.size as usize;
        leaf.hashes[..=n].copy_from_slice(&other.hashes[..=n]);
        leaf.entries[..n].copy_from_slice(&other.entries[..n]);
        leaf
    }

    /// HiGHS: HighsHashTree.h::InnerLeaf::insert_entry (lines 131-163).
    /// Returns `true` when inserted; the C++ value pointer is re-found by the caller.
    fn insert_entry(
        &mut self,
        full_hash: u64,
        hash_pos: i32,
        entry: HighsHashTableEntry<K, V>,
    ) -> bool {
        debug_assert!(self.size < Self::capacity());
        let hash: u16 = get_hash_chunks16(full_hash, hash_pos);
        let hash_chunk: u8 = get_first_chunk16(hash);

        let mut pos = self.occupation.num_set_until(hash_chunk);

        if self.occupation.test(hash_chunk) {
            // since the occupation flag is set we need to start searching from
            // pos-1 and can rely on a hash chunk with the same value existing for
            // the scan
            pos -= 1;
            while self.hashes[pos as usize] > hash as u64 {
                pos += 1;
            }

            if self.find_key(entry.key(), hash, &mut pos) {
                return false;
            }
        } else {
            self.occupation.set(hash_chunk);

            if pos < self.size {
                while self.hashes[pos as usize] > hash as u64 {
                    pos += 1;
                }
            }
        }

        if pos < self.size {
            let size = self.size;
            self.move_backward(pos, size);
        }
        self.entries[pos as usize] = entry;
        self.hashes[pos as usize] = hash as u64;
        self.size += 1;
        self.hashes[self.size as usize] = 0;

        true
    }

    /// The scan shared by both `find_entry` overloads: the slot of `key`, if present.
    #[inline]
    fn find_entry_pos(&self, full_hash: u64, hash_pos: i32, key: &K) -> Option<usize> {
        let hash: u16 = get_hash_chunks16(full_hash, hash_pos);
        let hash_chunk: u8 = get_first_chunk16(hash);
        if !self.occupation.test(hash_chunk) {
            return None;
        }

        let mut pos = self.occupation.num_set_until(hash_chunk) - 1;
        while self.hashes[pos as usize] > hash as u64 {
            pos += 1;
        }

        if self.find_key(key, hash, &mut pos) {
            return Some(pos as usize);
        }

        None
    }

    /// HiGHS: HighsHashTree.h::InnerLeaf::find_entry (lines 165-176)
    #[inline]
    fn find_entry(&self, full_hash: u64, hash_pos: i32, key: &K) -> Option<&V> {
        let pos = self.find_entry_pos(full_hash, hash_pos, key)?;
        Some(self.entries[pos].value())
    }

    /// HiGHS: HighsHashTree.h::InnerLeaf::find_entry (lines 165-176), non-const overload.
    #[inline]
    fn find_entry_mut(&mut self, full_hash: u64, hash_pos: i32, key: &K) -> Option<&mut V> {
        let pos = self.find_entry_pos(full_hash, hash_pos, key)?;
        Some(&mut self.entries[pos].value_)
    }

    /// HiGHS: HighsHashTree.h::InnerLeaf::erase_entry (lines 178-205)
    fn erase_entry(&mut self, full_hash: u64, hash_pos: i32, key: &K) -> bool {
        let hash: u16 = get_hash_chunks16(full_hash, hash_pos);
        let hash_chunk: u8 = get_first_chunk16(hash);
        if !self.occupation.test(hash_chunk) {
            return false;
        }

        let mut start_pos = self.occupation.num_set_until(hash_chunk) - 1;
        while get_first_chunk16(self.hashes[start_pos as usize] as u16) > hash_chunk {
            start_pos += 1;
        }

        let mut pos = start_pos;
        while self.hashes[pos as usize] > hash as u64 {
            pos += 1;
        }

        if !self.find_key(key, hash, &mut pos) {
            return false;
        }

        self.size -= 1;
        if pos < self.size {
            // std::move(entries[pos+1 .. size+1), entries.begin()+pos)
            self.entries
                .copy_within((pos + 1) as usize..(self.size + 1) as usize, pos as usize);
            // memmove(&hashes[pos], &hashes[pos+1], sizeof * (size - pos))
            self.hashes
                .copy_within((pos + 1) as usize..(self.size + 1) as usize, pos as usize);
            if get_first_chunk16(self.hashes[start_pos as usize] as u16) != hash_chunk {
                self.occupation.flip(hash_chunk);
            }
        } else if start_pos == pos {
            self.occupation.flip(hash_chunk);
        }

        self.hashes[self.size as usize] = 0;
        true
    }

    /// HiGHS: HighsHashTree.h::InnerLeaf::rehash (lines 207-255)
    fn rehash(&mut self, hash_pos: i32) {
        // function needs to possibly reorder elements by a different hash value
        // (comment kept from C++: chances are very high we are already ordered)
        self.occupation = HighsHashTree_Occupation::default();
        for i in 0..self.size as usize {
            let full = self.entries[i].key().highs_hash();
            self.hashes[i] = get_hash_chunks16(full, hash_pos) as u64;
            self.occupation
                .set(get_first_chunk16(self.hashes[i] as u16));
        }

        let mut i: i32 = 0;
        while i < self.size {
            let hash_chunk = get_first_chunk16(self.hashes[i as usize] as u16);
            let mut pos = self.occupation.num_set_until(hash_chunk) - 1;

            // if the position is after i the element definitely comes later, so we
            // swap it to that position and proceed without increasing i until
            // eventually an element appears that comes at position i or before
            if pos > i {
                self.hashes.swap(pos as usize, i as usize);
                self.entries.swap(pos as usize, i as usize);
                continue;
            }

            // the position is before or at i, now check where the exact location
            // should be for the ordering by hash so that the invariant is that all
            // elements up to i are properly sorted.
            while pos < i && self.hashes[pos as usize] >= self.hashes[i as usize] {
                pos += 1;
            }

            // if the final position is before i we need to move elements to
            // make space at that position, otherwise nothing needs to be done but
            // incrementing i increasing the sorted range by 1.
            if pos < i {
                let hash = self.hashes[i as usize];
                let entry = self.entries[i as usize];
                self.move_backward(pos, i);
                self.hashes[pos as usize] = hash;
                self.entries[pos as usize] = entry;
            }
            i += 1;
        }
    }

    /// HiGHS: HighsHashTree.h::InnerLeaf::move_backward (lines 257-264)
    #[inline]
    fn move_backward(&mut self, first: i32, last: i32) {
        // std::move_backward(entries[first..last) -> ..., entries.begin()+last+1); the
        // entries are trivially copyable, so this is the memmove C++ compiles it to
        self.entries
            .copy_within(first as usize..last as usize, (first + 1) as usize);
        // memmove(&hashes[first + 1], &hashes[first], sizeof * (last - first))
        self.hashes
            .copy_within(first as usize..last as usize, (first + 1) as usize);
    }

    /// HiGHS: HighsHashTree.h::InnerLeaf::find_key (lines 266-273)
    #[inline]
    fn find_key(&self, key: &K, hash: u16, pos: &mut i32) -> bool {
        // find key
        while *pos != self.size && self.hashes[*pos as usize] == hash as u64 {
            if key == self.entries[*pos as usize].key() {
                return true;
            }
            *pos += 1;
        }
        false
    }

    /// The live entries `entries[0..size)`.
    #[inline]
    fn live(&self) -> &[HighsHashTableEntry<K, V>] {
        &self.entries[..self.size as usize]
    }

    /// HiGHS 1.16 #3327 (h16_hash_common): HighsHashTree.h::HighsHashTree::pastChunkEnd
    /// (lines 507-510 of 1.16.0)
    #[inline]
    fn past_chunk_end(&self, i: i32, pos: i32) -> bool {
        i >= self.size || get_first_chunk16(self.hashes[i as usize] as u16) as i32 != pos
    }

    /// HiGHS 1.16 #3327 (h16_hash_common): HighsHashTree.h::HighsHashTree::advancePastCurrentHash
    /// (lines 512-518 of 1.16.0)
    #[inline]
    fn advance_past_current_hash(&self, i: &mut i32, pos: i32) {
        let common_hash = self.hashes[*i as usize] as u16;
        *i += 1;
        while !self.past_chunk_end(*i, pos) && self.hashes[*i as usize] == common_hash as u64 {
            *i += 1;
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::mergeIntoLeaf (lines 461-506), the target
    /// leaf being `self`.
    fn merge_from(&mut self, hash_pos: i32, node: HighsHashTreeNode<K, V>) {
        match node {
            HighsHashTreeNode::ListLeaf(merge_leaf) => {
                let HighsHashTree_ListLeaf { first, .. } = *merge_leaf;
                let HighsHashTree_ListNode { next, entry } = first;
                self.insert_entry(entry.key().highs_hash(), hash_pos, entry);
                let mut iter = next;
                while let Some(node) = iter {
                    let HighsHashTree_ListNode { next, entry } = *node;
                    self.insert_entry(entry.key().highs_hash(), hash_pos, entry);
                    iter = next;
                }
            }
            // C++ mergeIntoLeaf(leaf, mergeLeaf, hashPos) for the inner leaf classes
            other => with_inner_leaf!(other, merge_leaf => {
                for &entry in merge_leaf.live() {
                    self.insert_entry(entry.key().highs_hash(), hash_pos, entry);
                }
            }, _ => {}),
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::findCommonInLeaf, leaf-leaf overload
    /// (lines 508-558; 1.16.0 lines 520-570 for `h16_hash_common`).
    fn find_common_with_leaf<'t, const CAP2: usize, const HCAP2: usize>(
        &'t self,
        leaf2: &HighsHashTree_InnerLeaf<K, V, CAP2, HCAP2>,
        h16_hash_common: bool,
    ) -> Option<&'t HighsHashTableEntry<K, V>> {
        let leaf = self;
        let mut match_mask = leaf.occupation.occupation & leaf2.occupation.occupation;
        if match_mask == 0 {
            return None;
        }

        let mut offset1: i32 = -1;
        let mut offset2: i32 = -1;
        while match_mask != 0 {
            let pos = crate::util::highs_hash::log2i_u64(match_mask);
            match_mask ^= 1u64 << pos;

            let mut i = leaf.occupation.num_set_until(pos as u8) + offset1;
            while get_first_chunk16(leaf.hashes[i as usize] as u16) as i32 != pos {
                i += 1;
                offset1 += 1;
            }

            let mut j = leaf2.occupation.num_set_until(pos as u8) + offset2;
            while get_first_chunk16(leaf2.hashes[j as usize] as u16) as i32 != pos {
                j += 1;
                offset2 += 1;
            }

            if h16_hash_common {
                // HiGHS 1.16 #3327 (h16_hash_common): HighsHashTree.h lines 546-562 (1.16.0)
                loop {
                    if leaf.hashes[i as usize] > leaf2.hashes[j as usize] {
                        i += 1;
                    } else if leaf2.hashes[j as usize] > leaf.hashes[i as usize] {
                        j += 1;
                    } else {
                        let i_save = i;
                        let j_save = j;
                        leaf.advance_past_current_hash(&mut i, pos);
                        leaf2.advance_past_current_hash(&mut j, pos);
                        for ii in i_save..i {
                            for jj in j_save..j {
                                if leaf.entries[ii as usize].key()
                                    == leaf2.entries[jj as usize].key()
                                {
                                    return Some(&leaf.entries[ii as usize]);
                                }
                            }
                        }
                    }
                    if leaf.past_chunk_end(i, pos) || leaf2.past_chunk_end(j, pos) {
                        break;
                    }
                }
            } else {
                loop {
                    if leaf.hashes[i as usize] > leaf2.hashes[j as usize] {
                        i += 1;
                        if i == leaf.size
                            || get_first_chunk16(leaf.hashes[i as usize] as u16) as i32 != pos
                        {
                            break;
                        }
                    } else if leaf2.hashes[j as usize] > leaf.hashes[i as usize] {
                        j += 1;
                        if j == leaf2.size
                            || get_first_chunk16(leaf2.hashes[j as usize] as u16) as i32 != pos
                        {
                            break;
                        }
                    } else {
                        if leaf.entries[i as usize].key() == leaf2.entries[j as usize].key() {
                            return Some(&leaf.entries[i as usize]);
                        }

                        i += 1;
                        if i == leaf.size
                            || get_first_chunk16(leaf.hashes[i as usize] as u16) as i32 != pos
                        {
                            break;
                        }
                        j += 1;
                        if j == leaf2.size
                            || get_first_chunk16(leaf2.hashes[j as usize] as u16) as i32 != pos
                        {
                            break;
                        }
                    }
                }
            }
        }

        None
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::findCommonInLeaf, leaf-branch overload
    /// (lines 560-605; 1.16.0 lines 572-614 for `h16_hash_common`).
    fn find_common_with_branch<'t>(
        &'t self,
        branch: &HighsHashTree_BranchNode<K, V>,
        hash_pos: i32,
        h16_hash_common: bool,
    ) -> Option<&'t HighsHashTableEntry<K, V>> {
        let leaf = self;
        let mut match_mask = branch.occupation.occupation & leaf.occupation.occupation;

        let mut offset: i32 = -1;
        while match_mask != 0 {
            let pos = crate::util::highs_hash::log2i_u64(match_mask);
            match_mask ^= 1u64 << pos;

            let mut i = leaf.occupation.num_set_until(pos as u8) + offset;
            while get_first_chunk16(leaf.hashes[i as usize] as u16) as i32 != pos {
                i += 1;
                offset += 1;
            }

            let j = branch.occupation.num_set_until(pos as u8) - 1;

            loop {
                let key = leaf.entries[i as usize].key();
                if HighsHashTree::<K, V>::find_recurse(
                    &branch.children[j as usize],
                    key.highs_hash(),
                    hash_pos + 1,
                    key,
                )
                .is_some()
                {
                    return Some(&leaf.entries[i as usize]);
                }
                i += 1;
                // HiGHS 1.16 #3327 (h16_hash_common): HighsHashTree.h line 605 (1.16.0)
                let more = if h16_hash_common {
                    !leaf.past_chunk_end(i, pos)
                } else {
                    i < leaf.size && get_first_chunk16(leaf.hashes[i as usize] as u16) as i32 == pos
                };
                if !more {
                    break;
                }
            }
        }
        None
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::findCommonInLeaf (lines 508-605)
    /// HiGHS 1.16 #3327 (h16_hash_common): HighsHashTree.h::findCommonInLeaf (lines 507-614 of 1.16.0)
    fn find_common_in_leaf<'t>(
        &'t self,
        n2: &HighsHashTreeNode<K, V>,
        hash_pos: i32,
        h16_hash_common: bool,
    ) -> Option<&'t HighsHashTableEntry<K, V>> {
        match n2 {
            HighsHashTreeNode::Branch(branch) => {
                self.find_common_with_branch(branch, hash_pos, h16_hash_common)
            }
            other => with_inner_leaf!(other, leaf2 => {
                self.find_common_with_leaf(leaf2, h16_hash_common)
            }, _ => None),
        }
    }
}

/// HiGHS: HighsHashTree.h::HighsHashTree::get_hash_chunks16 (lines 55-57)
#[inline]
fn get_hash_chunks16(hash: u64, pos: i32) -> u16 {
    (hash >> (48 - pos * kBitsPerLevel)) as u16
}

/// HiGHS: HighsHashTree.h::HighsHashTree::get_first_chunk16 (lines 59-61)
#[inline]
fn get_first_chunk16(chunks: u16) -> u8 {
    (chunks >> (16 - kBitsPerLevel)) as u8
}

/// A new, empty boxed inner leaf of class `size_class` (1..=4) as a node
/// (C++ `new InnerLeaf<kSizeClass>` in the `switch (sizeClass)` blocks).
fn new_inner_leaf_node<K: HighsHashKey + Copy, V: Copy + Default>(
    size_class: i32,
) -> HighsHashTreeNode<K, V> {
    match size_class {
        1 => HighsHashTreeNode::InnerLeaf1(HighsHashTree_InnerLeaf1::new_boxed()),
        2 => HighsHashTreeNode::InnerLeaf2(HighsHashTree_InnerLeaf2::new_boxed()),
        3 => HighsHashTreeNode::InnerLeaf3(HighsHashTree_InnerLeaf3::new_boxed()),
        4 => HighsHashTreeNode::InnerLeaf4(HighsHashTree_InnerLeaf4::new_boxed()),
        _ => {
            // Unexpected result from 'entries_to_size_class'
            debug_assert!(false, "Unexpected result from 'entries_to_size_class'");
            HighsHashTreeNode::Empty
        }
    }
}

impl<K: HighsHashKey + Copy, V: Copy + Default> HighsHashTree<K, V> {
    /// HiGHS: HighsHashTree.h::HighsHashTree::compute_hash (lines 46-48)
    pub fn compute_hash(key: &K) -> u64 {
        key.highs_hash()
    }
    /// HiGHS: HighsHashTree.h::HighsHashTree::get_hash_chunk (lines 50-53)
    pub fn get_hash_chunk(hash: u64, pos: i32) -> u8 {
        ((hash >> (64 - kBitsPerLevel - pos * kBitsPerLevel)) & ((kBranchFactor - 1) as u64)) as u8
    }
    /// HiGHS: HighsHashTree.h::HighsHashTree::get_hash_chunks16 (lines 55-57)
    pub fn get_hash_chunks16(hash: u64, pos: i32) -> u16 {
        (hash >> (48 - pos * kBitsPerLevel)) as u16
    }
    /// HiGHS: HighsHashTree.h::HighsHashTree::get_first_chunk16 (lines 59-61)
    pub fn get_first_chunk16(chunks: u16) -> u8 {
        (chunks >> (16 - kBitsPerLevel)) as u8
    }
    /// HiGHS: HighsHashTree.h::HighsHashTree::set_hash_chunk (lines 63-67)
    pub fn set_hash_chunk(hash: &mut u64, chunk: u64, chunk_pos: i32) {
        let shift_amount = 60 - kBitsPerLevel - chunk_pos * kBitsPerLevel;
        let chunk = chunk ^ ((*hash >> shift_amount) & ((kBranchFactor - 1) as u64));
        *hash ^= chunk << shift_amount;
    }
    /// HiGHS: HighsHashTree.h::HighsHashTree::entries_to_size_class (lines 93-97)
    pub fn entries_to_size_class(num_entries: u32) -> i32 {
        let per_class: u32 = ((kLeafBurstThreshold - kMinLeafSize) / 3) as u32;
        1 + ((num_entries + (per_class - kMinLeafSize as u32 - 1)) / per_class) as i32
    }
    /// HiGHS: HighsHashTree.h::HighsHashTree::addChildToBranchNode (lines 431-459)
    /// C++ also sets the new child slot to nullptr in the caller; here the inserted slot is `Empty`.
    pub fn add_child_to_branch_node(
        branch: &mut HighsHashTree_BranchNode<K, V>,
        hash_value: u8,
        location: i32,
    ) {
        let right_childs = branch.occupation.num_set_after(hash_value);
        debug_assert!(right_childs + location == branch.occupation.num_set());
        branch
            .children
            .insert(location as usize, HighsHashTreeNode::Empty);
    }
    /// HiGHS: HighsHashTree.h::HighsHashTree::removeChildFromBranchNode (lines 607-705)
    ///
    /// Takes the branch by value: C++ frees or shrinks it and returns the
    /// replacement node. The removed child is still in `children` on entry (the
    /// caller has already cleared its occupation bit), exactly as in C++.
    #[allow(clippy::boxed_local)]
    pub fn remove_child_from_branch_node(
        branch: Box<HighsHashTree_BranchNode<K, V>>,
        location: i32,
        hash: u64,
        hash_pos: i32,
    ) -> HighsHashTreeNode<K, V> {
        let mut branch = branch;
        let new_num_child = branch.occupation.num_set();

        // first check if we might be able to merge all children into one new leaf
        // based on the node numbers and assuming all of them might be in the
        // smallest size class
        if new_num_child * inner_leaf_capacity(1) <= kLeafBurstThreshold {
            // since we have a good chance of merging we now check the actual size
            // classes to see if that yields a number of entries at most the burst
            // threshold
            let mut child_entries: i32 = 0;
            for i in 0..=(new_num_child as usize) {
                child_entries += branch.children[i].num_entries_estimate();
                if child_entries > kLeafBurstThreshold {
                    break;
                }
            }

            if child_entries < kLeafBurstThreshold {
                // recompute the number of entries from the actual children
                child_entries = 0;
                for i in 0..=(new_num_child as usize) {
                    child_entries += branch.children[i].num_entries();
                }

                // check again if we exceed due to the extremely unlikely case
                // of having less than 5 list nodes with together more than 30 entries
                // as list nodes are only created in the last depth level
                if child_entries < kLeafBurstThreshold {
                    let size_class = Self::entries_to_size_class(child_entries as u32);
                    // an unexpected size class yields `Empty` (newNode stays kEmpty, as in C++)
                    let mut new_node = new_inner_leaf_node(size_class);
                    with_inner_leaf!(&mut new_node, new_leaf => {
                        for i in 0..=(new_num_child as usize) {
                            let child = std::mem::take(&mut branch.children[i]);
                            new_leaf.merge_from(hash_pos, child);
                        }
                    }, _ => {});
                    // destroyBranchingNode(branch) is the drop of `branch` here
                    return new_node;
                }
            }
        }

        // C++ memmove(childPtr(location), childPtr(location + 1), ...): drop the empty
        // child slot. The allocation size class change of C++ has no Rust counterpart.
        let _ = hash;
        branch.children.remove(location as usize);
        HighsHashTreeNode::Branch(branch)
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::insert_into_leaf (lines 707-722)
    /// `insert_node` is the leaf node (size class 1..3) itself; a full leaf of class
    /// `kSizeClass` is moved into a new leaf of class `kSizeClass + 1`.
    pub fn insert_into_leaf(
        insert_node: &mut HighsHashTreeNode<K, V>,
        hash: u64,
        hash_pos: i32,
        entry: HighsHashTableEntry<K, V>,
    ) -> bool {
        macro_rules! grow {
            ($leaf:ident, $Next:ident, $NextLeaf:ident) => {{
                if $leaf.size == $leaf.entries.len() as i32 {
                    if $leaf.find_entry(hash, hash_pos, entry.key()).is_some() {
                        return false;
                    }
                    let mut new_leaf = $NextLeaf::from_moved_leaf(&**$leaf);
                    let inserted = new_leaf.insert_entry(hash, hash_pos, entry);
                    *insert_node = HighsHashTreeNode::$Next(new_leaf);
                    return inserted;
                }
                $leaf.insert_entry(hash, hash_pos, entry)
            }};
        }
        match insert_node {
            HighsHashTreeNode::InnerLeaf1(leaf) => {
                grow!(leaf, InnerLeaf2, HighsHashTree_InnerLeaf2)
            }
            HighsHashTreeNode::InnerLeaf2(leaf) => {
                grow!(leaf, InnerLeaf3, HighsHashTree_InnerLeaf3)
            }
            HighsHashTreeNode::InnerLeaf3(leaf) => {
                grow!(leaf, InnerLeaf4, HighsHashTree_InnerLeaf4)
            }
            _ => unreachable!("insert_into_leaf called on a non-leaf node"),
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::insert_recurse (lines 724-958)
    /// Returns `inserted`; `insert_or_get` re-finds the value (C++ returns `pair<ValueType*, bool>`).
    /// The C++ tail call `return insert_recurse(insertNode, ...)` is the `loop` below.
    pub fn insert_recurse(
        insert_node: &mut HighsHashTreeNode<K, V>,
        hash: u64,
        hash_pos: i32,
        entry: HighsHashTableEntry<K, V>,
    ) -> bool {
        let mut hash_pos = hash_pos;
        let mut cur: &mut HighsHashTreeNode<K, V> = insert_node;
        loop {
            match cur {
                HighsHashTreeNode::Empty => {
                    if hash_pos == kMaxDepth {
                        *cur = HighsHashTreeNode::ListLeaf(Box::new(HighsHashTree_ListLeaf {
                            first: HighsHashTree_ListNode { next: None, entry },
                            count: 1,
                        }));
                        return true;
                    } else {
                        let mut leaf = HighsHashTree_InnerLeaf1::new_boxed();
                        let inserted = leaf.insert_entry(hash, hash_pos, entry);
                        *cur = HighsHashTreeNode::InnerLeaf1(leaf);
                        return inserted;
                    }
                }
                HighsHashTreeNode::ListLeaf(leaf) => {
                    let HighsHashTree_ListLeaf { first, count } = &mut **leaf;
                    let mut iter: &mut HighsHashTree_ListNode<K, V> = first;
                    loop {
                        // check for existing key
                        if iter.entry.key() == entry.key() {
                            return false;
                        }

                        if iter.next.is_none() {
                            // reached the end of the list and key is not duplicate, so insert
                            iter.next =
                                Some(Box::new(HighsHashTree_ListNode { next: None, entry }));
                            *count += 1;
                            return true;
                        }
                        iter = iter.next.as_deref_mut().unwrap();
                    }
                }
                HighsHashTreeNode::InnerLeaf1(_)
                | HighsHashTreeNode::InnerLeaf2(_)
                | HighsHashTreeNode::InnerLeaf3(_) => {
                    return Self::insert_into_leaf(&mut *cur, hash, hash_pos, entry);
                }
                HighsHashTreeNode::InnerLeaf4(leaf) => {
                    // kInnerLeafSizeClass4
                    if leaf.size < HighsHashTree_InnerLeaf4::<K, V>::capacity() {
                        return leaf.insert_entry(hash, hash_pos, entry);
                    }

                    if leaf.find_entry(hash, hash_pos, entry.key()).is_some() {
                        return false;
                    }
                    let mut occupation = leaf.occupation;

                    let hash_chunk = Self::get_hash_chunk(hash, hash_pos);
                    occupation.set(hash_chunk);

                    let branch_size = occupation.num_set();

                    // createBranchingNode(branchSize): C++ leaves the children uninitialised,
                    // they are set below (`Empty` here).
                    let mut branch = Box::new(HighsHashTree_BranchNode {
                        occupation: HighsHashTree_Occupation::default(),
                        children: (0..branch_size).map(|_| HighsHashTreeNode::Empty).collect(),
                    });
                    branch.occupation = occupation;

                    // C++ `*insertNode = branch` happens before the old leaf is freed;
                    // here the old leaf is taken out and `*cur` is written at each exit.
                    let mut old = match std::mem::take(&mut *cur) {
                        HighsHashTreeNode::InnerLeaf4(old) => old,
                        _ => unreachable!(),
                    };

                    if hash_pos + 1 == kMaxDepth {
                        // children are all null (Empty) at this point
                        for i in 0..old.size as usize {
                            let pos = occupation
                                .num_set_until(Self::get_first_chunk16(old.hashes[i] as u16))
                                - 1;
                            let e = old.entries[i];
                            if branch.children[pos as usize].get_type() == kEmpty {
                                branch.children[pos as usize] =
                                    HighsHashTreeNode::ListLeaf(Box::new(HighsHashTree_ListLeaf {
                                        first: HighsHashTree_ListNode {
                                            next: None,
                                            entry: e,
                                        },
                                        count: 1,
                                    }));
                            } else {
                                match &mut branch.children[pos as usize] {
                                    HighsHashTreeNode::ListLeaf(list) => {
                                        // newNode = new ListNode(move(first)); first.next = newNode;
                                        // first.entry = move(entry); ++count
                                        let old_first = std::mem::replace(
                                            &mut list.first,
                                            HighsHashTree_ListNode {
                                                next: None,
                                                entry: e,
                                            },
                                        );
                                        list.first.next = Some(Box::new(old_first));
                                        list.count += 1;
                                    }
                                    _ => unreachable!(),
                                }
                            }
                        }

                        drop(old);

                        let pos =
                            occupation.num_set_until(Self::get_hash_chunk(hash, hash_pos)) - 1;
                        if branch.children[pos as usize].get_type() == kEmpty {
                            branch.children[pos as usize] =
                                HighsHashTreeNode::ListLeaf(Box::new(HighsHashTree_ListLeaf {
                                    first: HighsHashTree_ListNode { next: None, entry },
                                    count: 1,
                                }));
                        } else {
                            match &mut branch.children[pos as usize] {
                                HighsHashTreeNode::ListLeaf(list) => {
                                    let old_first = std::mem::replace(
                                        &mut list.first,
                                        HighsHashTree_ListNode { next: None, entry },
                                    );
                                    list.first.next = Some(Box::new(old_first));
                                    list.count += 1;
                                }
                                _ => unreachable!(),
                            }
                        }
                        *cur = HighsHashTreeNode::Branch(branch);
                        return true;
                    }

                    if branch_size > 1 {
                        // maxsize in one bucket = number of items - (num buckets-1)
                        // since each bucket has at least 1 item the largest one can only
                        // have all remaining ones After adding the item: If it does not
                        // collide
                        let max_entries_per_leaf = 2 + old.size - branch_size;

                        if max_entries_per_leaf <= HighsHashTree_InnerLeaf1::<K, V>::capacity() {
                            // all items can go into the smallest leaf size
                            for i in 0..branch_size as usize {
                                branch.children[i] = HighsHashTreeNode::InnerLeaf1(
                                    HighsHashTree_InnerLeaf1::new_boxed(),
                                );
                            }

                            for i in 0..old.size as usize {
                                let pos = occupation
                                    .num_set_until(Self::get_first_chunk16(old.hashes[i] as u16))
                                    - 1;
                                let e = old.entries[i];
                                let h = Self::compute_hash(e.key());
                                match &mut branch.children[pos as usize] {
                                    HighsHashTreeNode::InnerLeaf1(cl) => {
                                        cl.insert_entry(h, hash_pos + 1, e);
                                    }
                                    _ => unreachable!(),
                                }
                            }

                            drop(old);

                            let pos =
                                occupation.num_set_until(Self::get_hash_chunk(hash, hash_pos)) - 1;
                            let inserted = match &mut branch.children[pos as usize] {
                                HighsHashTreeNode::InnerLeaf1(cl) => {
                                    cl.insert_entry(hash, hash_pos + 1, entry)
                                }
                                _ => unreachable!(),
                            };
                            *cur = HighsHashTreeNode::Branch(branch);
                            return inserted;
                        } else {
                            // there are many collisions, determine the exact sizes first
                            let mut sizes = [0u8; 55];
                            sizes[(occupation.num_set_until(hash_chunk) - 1) as usize] += 1;
                            for i in 0..old.size as usize {
                                let pos = occupation
                                    .num_set_until(Self::get_first_chunk16(old.hashes[i] as u16))
                                    - 1;
                                sizes[pos as usize] += 1;
                            }

                            for i in 0..branch_size as usize {
                                let size_class = Self::entries_to_size_class(sizes[i] as u32);
                                branch.children[i] = new_inner_leaf_node(size_class);
                            }

                            for i in 0..old.size as usize {
                                let pos = occupation
                                    .num_set_until(Self::get_first_chunk16(old.hashes[i] as u16))
                                    - 1;
                                let e = old.entries[i];
                                let h = Self::compute_hash(e.key());
                                // C++ default: break (nothing to insert into)
                                with_inner_leaf!(&mut branch.children[pos as usize], cl => {
                                    cl.insert_entry(h, hash_pos + 1, e);
                                }, _ => {});
                            }
                        }

                        drop(old);

                        let pos = occupation.num_set_until(hash_chunk) - 1;
                        hash_pos += 1;
                        *cur = HighsHashTreeNode::Branch(branch);
                        let HighsHashTreeNode::Branch(b) = cur else {
                            unreachable!()
                        };
                        cur = &mut b.children[pos as usize];
                    } else {
                        // extremely unlikely that the new branch node only gets one
                        // child in that case create it and defer the insertion into
                        // the next depth
                        hash_pos += 1;
                        old.rehash(hash_pos);
                        branch.children[0] = HighsHashTreeNode::InnerLeaf4(old);
                        *cur = HighsHashTreeNode::Branch(branch);
                        let HighsHashTreeNode::Branch(b) = cur else {
                            unreachable!()
                        };
                        cur = &mut b.children[0];
                    }
                }
                HighsHashTreeNode::Branch(branch) => {
                    let chunk = Self::get_hash_chunk(hash, hash_pos);
                    let mut location = branch.occupation.num_set_until(chunk);

                    if branch.occupation.test(chunk) {
                        location -= 1;
                    } else {
                        Self::add_child_to_branch_node(branch, chunk, location);
                        branch.occupation.set(chunk);
                    }

                    hash_pos += 1;
                    cur = &mut branch.children[location as usize];
                }
            }
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::erase_recurse (lines 960-1066)
    pub fn erase_recurse(
        erase_node: &mut HighsHashTreeNode<K, V>,
        hash: u64,
        hash_pos: i32,
        key: &K,
    ) {
        match erase_node {
            HighsHashTreeNode::Empty => {}
            HighsHashTreeNode::ListLeaf(leaf) => {
                let HighsHashTree_ListLeaf { first, count } = &mut **leaf;

                // check for existing key
                let mut iter: &mut HighsHashTree_ListNode<K, V> = first;
                loop {
                    if iter.entry.key() == key {
                        // key found, decrease count
                        *count -= 1;
                        // if we have a next node after replace the entry in iter by
                        // moving that node into it (C++ `*iter = std::move(*next)`)
                        if let Some(next) = iter.next.take() {
                            let HighsHashTree_ListNode {
                                next: next_next,
                                entry,
                            } = *next;
                            iter.entry = entry;
                            iter.next = next_next;
                        }
                        // HiGHS quirk (verified against HighsHashTree.h:966-995): when the
                        // matching node is the last one and not `first` (no `next`), C++ only
                        // decrements `count`; the predecessor still links it, so the key stays
                        // findable (and blocks a re-insert). Reproduced faithfully. List leaves
                        // exist only at kMaxDepth (54 equal hash bits), so this needs a collision.
                        break;
                    }

                    match iter.next.as_deref_mut() {
                        Some(n) => iter = n,
                        None => break,
                    }
                }

                let remaining = *count;
                if remaining == 0 {
                    *erase_node = HighsHashTreeNode::Empty;
                }
            }
            HighsHashTreeNode::InnerLeaf1(leaf) => {
                if leaf.erase_entry(hash, hash_pos, key) && leaf.size == 0 {
                    *erase_node = HighsHashTreeNode::Empty;
                }
            }
            // a leaf of class kSizeClass > 1 that drops to the capacity of class
            // kSizeClass - 1 moves into a leaf of that class
            HighsHashTreeNode::InnerLeaf2(leaf) => {
                if leaf.erase_entry(hash, hash_pos, key)
                    && leaf.size == HighsHashTree_InnerLeaf1::<K, V>::capacity()
                {
                    *erase_node = HighsHashTreeNode::InnerLeaf1(
                        HighsHashTree_InnerLeaf1::from_moved_leaf(&**leaf),
                    );
                }
            }
            HighsHashTreeNode::InnerLeaf3(leaf) => {
                if leaf.erase_entry(hash, hash_pos, key)
                    && leaf.size == HighsHashTree_InnerLeaf2::<K, V>::capacity()
                {
                    *erase_node = HighsHashTreeNode::InnerLeaf2(
                        HighsHashTree_InnerLeaf2::from_moved_leaf(&**leaf),
                    );
                }
            }
            HighsHashTreeNode::InnerLeaf4(leaf) => {
                if leaf.erase_entry(hash, hash_pos, key)
                    && leaf.size == HighsHashTree_InnerLeaf3::<K, V>::capacity()
                {
                    *erase_node = HighsHashTreeNode::InnerLeaf3(
                        HighsHashTree_InnerLeaf3::from_moved_leaf(&**leaf),
                    );
                }
            }
            HighsHashTreeNode::Branch(branch) => {
                let chunk = Self::get_hash_chunk(hash, hash_pos);
                if !branch.occupation.test(chunk) {
                    return;
                }

                let location = branch.occupation.num_set_until(chunk) - 1;
                Self::erase_recurse(
                    &mut branch.children[location as usize],
                    hash,
                    hash_pos + 1,
                    key,
                );

                if branch.children[location as usize].get_type() != kEmpty {
                    return;
                }

                branch.occupation.flip(chunk);

                let old = match std::mem::take(&mut *erase_node) {
                    HighsHashTreeNode::Branch(b) => b,
                    _ => unreachable!(),
                };
                *erase_node = Self::remove_child_from_branch_node(old, location, hash, hash_pos);
            }
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::find_recurse (lines 1068-1113)
    pub fn find_recurse<'t>(
        node: &'t HighsHashTreeNode<K, V>,
        hash: u64,
        hash_pos: i32,
        key: &K,
    ) -> Option<&'t V> {
        let mut node = node;
        let mut hash_pos = hash_pos;
        loop {
            match node {
                HighsHashTreeNode::Empty => return None,
                HighsHashTreeNode::ListLeaf(leaf) => {
                    let mut iter: &'t HighsHashTree_ListNode<K, V> = &leaf.first;
                    loop {
                        if iter.entry.key() == key {
                            return Some(iter.entry.value());
                        }
                        iter = iter.next.as_deref()?;
                    }
                }
                HighsHashTreeNode::InnerLeaf1(leaf) => return leaf.find_entry(hash, hash_pos, key),
                HighsHashTreeNode::InnerLeaf2(leaf) => return leaf.find_entry(hash, hash_pos, key),
                HighsHashTreeNode::InnerLeaf3(leaf) => return leaf.find_entry(hash, hash_pos, key),
                HighsHashTreeNode::InnerLeaf4(leaf) => return leaf.find_entry(hash, hash_pos, key),
                HighsHashTreeNode::Branch(branch) => {
                    let chunk = Self::get_hash_chunk(hash, hash_pos);
                    if !branch.occupation.test(chunk) {
                        return None;
                    }
                    let location = branch.occupation.num_set_until(chunk) - 1;
                    node = &branch.children[location as usize];
                    hash_pos += 1;
                }
            }
        }
    }

    /// Private mirror of `find_recurse` returning a mutable value (for `find_mut` and
    /// `insert_or_get`). Same walk as `find_recurse`.
    fn find_recurse_mut<'t>(
        node: &'t mut HighsHashTreeNode<K, V>,
        hash: u64,
        hash_pos: i32,
        key: &K,
    ) -> Option<&'t mut V> {
        let mut node = node;
        let mut hash_pos = hash_pos;
        loop {
            match node {
                HighsHashTreeNode::Empty => return None,
                HighsHashTreeNode::ListLeaf(leaf) => {
                    // locate the matching node first, then re-walk it mutably
                    let mut depth = 0usize;
                    let mut found = false;
                    {
                        let mut iter: &HighsHashTree_ListNode<K, V> = &leaf.first;
                        loop {
                            if iter.entry.key() == key {
                                found = true;
                                break;
                            }
                            match iter.next.as_deref() {
                                Some(n) => {
                                    iter = n;
                                    depth += 1;
                                }
                                None => break,
                            }
                        }
                    }
                    if !found {
                        return None;
                    }
                    let mut iter: &'t mut HighsHashTree_ListNode<K, V> = &mut leaf.first;
                    for _ in 0..depth {
                        iter = iter.next.as_deref_mut().unwrap();
                    }
                    return Some(&mut iter.entry.value_);
                }
                HighsHashTreeNode::InnerLeaf1(leaf) => {
                    return leaf.find_entry_mut(hash, hash_pos, key);
                }
                HighsHashTreeNode::InnerLeaf2(leaf) => {
                    return leaf.find_entry_mut(hash, hash_pos, key);
                }
                HighsHashTreeNode::InnerLeaf3(leaf) => {
                    return leaf.find_entry_mut(hash, hash_pos, key);
                }
                HighsHashTreeNode::InnerLeaf4(leaf) => {
                    return leaf.find_entry_mut(hash, hash_pos, key);
                }
                HighsHashTreeNode::Branch(branch) => {
                    let chunk = Self::get_hash_chunk(hash, hash_pos);
                    if !branch.occupation.test(chunk) {
                        return None;
                    }
                    let location = branch.occupation.num_set_until(chunk) - 1;
                    node = &mut branch.children[location as usize];
                    hash_pos += 1;
                }
            }
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::find_common_recurse (lines 1115-1174)
    /// HiGHS 1.16 #3327 (h16_hash_common): new argument `h16_hash_common` selects the 1.16.0 findCommonInLeaf.
    pub fn find_common_recurse<'t>(
        n1: &'t HighsHashTreeNode<K, V>,
        n2: &'t HighsHashTreeNode<K, V>,
        hash_pos: i32,
        h16_hash_common: bool,
    ) -> Option<&'t HighsHashTableEntry<K, V>> {
        let mut n1 = n1;
        let mut n2 = n2;
        if n1.get_type() > n2.get_type() {
            std::mem::swap(&mut n1, &mut n2);
        }

        match n1 {
            HighsHashTreeNode::Empty => None,
            HighsHashTreeNode::ListLeaf(leaf) => {
                let mut iter: &'t HighsHashTree_ListNode<K, V> = &leaf.first;
                loop {
                    if Self::find_recurse(
                        n2,
                        Self::compute_hash(iter.entry.key()),
                        hash_pos,
                        iter.entry.key(),
                    )
                    .is_some()
                    {
                        return Some(&iter.entry);
                    }
                    iter = iter.next.as_deref()?;
                }
            }
            HighsHashTreeNode::InnerLeaf1(leaf1) => {
                leaf1.find_common_in_leaf(n2, hash_pos, h16_hash_common)
            }
            HighsHashTreeNode::InnerLeaf2(leaf1) => {
                leaf1.find_common_in_leaf(n2, hash_pos, h16_hash_common)
            }
            HighsHashTreeNode::InnerLeaf3(leaf1) => {
                leaf1.find_common_in_leaf(n2, hash_pos, h16_hash_common)
            }
            HighsHashTreeNode::InnerLeaf4(leaf1) => {
                leaf1.find_common_in_leaf(n2, hash_pos, h16_hash_common)
            }
            HighsHashTreeNode::Branch(branch1) => {
                let branch2 = match n2 {
                    HighsHashTreeNode::Branch(b) => b,
                    _ => unreachable!("Unexpected type in hash tree"),
                };

                let mut match_mask = branch1.occupation.occupation & branch2.occupation.occupation;

                while match_mask != 0 {
                    let pos = crate::util::highs_hash::log2i_u64(match_mask);
                    match_mask ^= 1u64 << pos;

                    let location1 = branch1.occupation.num_set_until(pos as u8) - 1;
                    let location2 = branch2.occupation.num_set_until(pos as u8) - 1;

                    let found = Self::find_common_recurse(
                        &branch1.children[location1 as usize],
                        &branch2.children[location2 as usize],
                        hash_pos + 1,
                        h16_hash_common,
                    );
                    if found.is_some() {
                        return found;
                    }
                }

                None
            }
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::for_each_recurse (lines 1266-1312; void `f`)
    pub fn for_each_recurse<F: FnMut(&K, &V) + ?Sized>(node: &HighsHashTreeNode<K, V>, f: &mut F) {
        match node {
            HighsHashTreeNode::Empty => {}
            HighsHashTreeNode::ListLeaf(leaf) => {
                let mut iter: &HighsHashTree_ListNode<K, V> = &leaf.first;
                loop {
                    f(iter.entry.key(), iter.entry.value());
                    match iter.next.as_deref() {
                        Some(n) => iter = n,
                        None => break,
                    }
                }
            }
            HighsHashTreeNode::Branch(branch) => {
                let size = branch.occupation.num_set();
                for i in 0..size as usize {
                    Self::for_each_recurse(&branch.children[i], &mut *f);
                }
            }
            other => with_inner_leaf!(other, leaf => {
                for e in leaf.live() {
                    f(e.key(), e.value());
                }
            }, _ => unreachable!()),
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::for_each_recurse (lines 1266-1312; void `f`),
    /// the form whose C++ lambda takes `V&`: `for_each` is `const` in C++, but it walks raw
    /// `NodePtr`s, so a lambda such as `[&](HighsInt, VarBound& vbd) { vbd.coef /= scale; }`
    /// (HighsImplications.h:142-150, HighsImplications.cpp:760-790) mutates the values.
    pub fn for_each_recurse_mut<F: FnMut(&K, &mut V) + ?Sized>(
        node: &mut HighsHashTreeNode<K, V>,
        f: &mut F,
    ) {
        match node {
            HighsHashTreeNode::Empty => {}
            HighsHashTreeNode::ListLeaf(leaf) => {
                let mut iter: &mut HighsHashTree_ListNode<K, V> = &mut leaf.first;
                loop {
                    let HighsHashTableEntry { key_, value_ } = &mut iter.entry;
                    f(key_, value_);
                    match iter.next.as_deref_mut() {
                        Some(n) => iter = n,
                        None => break,
                    }
                }
            }
            HighsHashTreeNode::Branch(branch) => {
                let size = branch.occupation.num_set();
                for i in 0..size as usize {
                    Self::for_each_recurse_mut(&mut branch.children[i], &mut *f);
                }
            }
            other => with_inner_leaf!(other, leaf => {
                let size = leaf.size as usize;
                for e in &mut leaf.entries[..size] {
                    let HighsHashTableEntry { key_, value_ } = e;
                    f(key_, value_);
                }
            }, _ => unreachable!()),
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::for_each_recurse (lines 1314-1378; `f` returns a
    /// truthy value that stops the iteration)
    pub fn for_each_recurse_until<F: FnMut(&K, &V) -> bool + ?Sized>(
        node: &HighsHashTreeNode<K, V>,
        f: &mut F,
    ) -> bool {
        match node {
            HighsHashTreeNode::Empty => false,
            HighsHashTreeNode::ListLeaf(leaf) => {
                let mut iter: &HighsHashTree_ListNode<K, V> = &leaf.first;
                loop {
                    if f(iter.entry.key(), iter.entry.value()) {
                        return true;
                    }
                    match iter.next.as_deref() {
                        Some(n) => iter = n,
                        None => return false,
                    }
                }
            }
            HighsHashTreeNode::Branch(branch) => {
                let size = branch.occupation.num_set();
                for i in 0..size as usize {
                    if Self::for_each_recurse_until(&branch.children[i], &mut *f) {
                        return true;
                    }
                }
                false
            }
            other => with_inner_leaf!(other, leaf => {
                for e in leaf.live() {
                    if f(e.key(), e.value()) {
                        return true;
                    }
                }
                false
            }, _ => unreachable!()),
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::insert (lines 1381-1386)
    /// C++: `template <typename... Args> bool insert(Args&&... args)` (builds an `Entry`)
    pub fn insert(&mut self, entry: HighsHashTableEntry<K, V>) -> bool {
        let hash = Self::compute_hash(entry.key());
        Self::insert_recurse(&mut self.root, hash, 0, entry)
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::insert_or_get (lines 1388-1393)
    pub fn insert_or_get(&mut self, entry: HighsHashTableEntry<K, V>) -> (&mut V, bool) {
        let hash = Self::compute_hash(entry.key());
        let key = *entry.key();
        let inserted = Self::insert_recurse(&mut self.root, hash, 0, entry);
        // the C++ value pointer is the entry just placed (or the existing one): re-find it
        match Self::find_recurse_mut(&mut self.root, hash, 0, &key) {
            Some(value) => (value, inserted),
            None => unreachable!("insert_or_get: key missing right after insert_recurse"),
        }
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::erase (lines 1395-1399)
    pub fn erase(&mut self, key: &K) {
        let hash = Self::compute_hash(key);
        Self::erase_recurse(&mut self.root, hash, 0, key);
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::contains (lines 1401-1404)
    pub fn contains(&self, key: &K) -> bool {
        let hash = Self::compute_hash(key);
        Self::find_recurse(&self.root, hash, 0, key).is_some()
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::find (lines 1406-1410)
    pub fn find(&self, key: &K) -> Option<&V> {
        let hash = Self::compute_hash(key);
        Self::find_recurse(&self.root, hash, 0, key)
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::find (lines 1412-1416, non-const)
    pub fn find_mut(&mut self, key: &K) -> Option<&mut V> {
        let hash = Self::compute_hash(key);
        Self::find_recurse_mut(&mut self.root, hash, 0, key)
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::find_common (lines 1418-1421)
    /// HiGHS 1.16 #3327 (h16_hash_common): new argument `h16_hash_common`, passed down to find_common_recurse.
    pub fn find_common<'a>(
        &'a self,
        other: &'a HighsHashTree<K, V>,
        h16_hash_common: bool,
    ) -> Option<&'a HighsHashTableEntry<K, V>> {
        Self::find_common_recurse(&self.root, &other.root, 0, h16_hash_common)
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::empty (line 1423)
    pub fn empty(&self) -> bool {
        self.root.get_type() == kEmpty
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::clear (lines 1425-1428)
    /// C++ `destroy_recurse(root); root = nullptr;`: dropping the old root frees it.
    pub fn clear(&mut self) {
        self.root = HighsHashTreeNode::Empty;
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::for_each (lines 1430-1435), `f` returning void.
    /// (`HighsHashTableEntry::forward(f)` calls `f(key)` for `V = ()` and `f(key, value)` otherwise:
    /// a key-only C++ lambda ignores the second argument here.)
    pub fn for_each<F: FnMut(&K, &V) + ?Sized>(&self, f: &mut F) {
        Self::for_each_recurse(&self.root, f);
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::for_each (lines 1430-1435), `f` taking `V&`
    /// (see [`Self::for_each_recurse_mut`]). Same visiting order as [`Self::for_each`].
    pub fn for_each_mut<F: FnMut(&K, &mut V) + ?Sized>(&mut self, f: &mut F) {
        Self::for_each_recurse_mut(&mut self.root, f);
    }

    /// HiGHS: HighsHashTree.h::HighsHashTree::for_each (lines 1430-1435), `f` returning `bool`:
    /// stops at the first `true` and returns it.
    pub fn for_each_until<F: FnMut(&K, &V) -> bool + ?Sized>(&self, f: &mut F) -> bool {
        Self::for_each_recurse_until(&self.root, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Randomized differential check against `HashMap` (shape-independent semantics:
    /// membership, values, counts and iteration totals).
    #[test]
    fn matches_hashmap_under_random_ops() {
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut tree: HighsHashTree<HighsInt, HighsInt> = HighsHashTree::default();
        let mut reference: HashMap<HighsInt, HighsInt> = HashMap::new();
        for step in 0..20_000u32 {
            let key = (next() % 3_000) as HighsInt;
            match next() % 3 {
                0 => {
                    let inserted = tree.insert(HighsHashTableEntry::new(key, step as HighsInt));
                    let expected = !reference.contains_key(&key);
                    assert_eq!(inserted, expected, "insert {key} at step {step}");
                    reference.entry(key).or_insert(step as HighsInt);
                }
                1 => {
                    tree.erase(&key);
                    reference.remove(&key);
                }
                _ => {
                    assert_eq!(
                        tree.find(&key).copied(),
                        reference.get(&key).copied(),
                        "find {key}"
                    );
                    assert_eq!(tree.contains(&key), reference.contains_key(&key));
                }
            }
        }
        let mut count = 0usize;
        tree.for_each(&mut |k, v| {
            assert_eq!(reference.get(k), Some(v));
            count += 1;
        });
        assert_eq!(count, reference.len());
        for (k, v) in &reference {
            assert_eq!(tree.find(k), Some(v));
        }
    }

    /// FNV-1a over the little-endian bytes, as the C++ generator computes it.
    struct Fnv(u64);
    impl Fnv {
        fn new() -> Self {
            Fnv(0xcbf2_9ce4_8422_2325)
        }
        fn add(&mut self, v: u64) {
            for i in 0..8 {
                self.0 ^= (v >> (8 * i)) & 0xff;
                self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
    }

    fn fnv_of(t: &HighsHashTree<HighsInt, HighsInt>) -> (i32, u64) {
        let mut f = Fnv::new();
        let mut n = 0;
        t.for_each(&mut |k, v| {
            f.add(u64::from(*k as u32));
            f.add(u64::from(*v as u32));
            n += 1;
        });
        (n, f.0)
    }

    /// Differential test against the real C++ `HighsHashTree` (generator and output in
    /// `tests/fixtures/mip/unit/hash_tree.{cpp,out}`). `for_each` order is the tree shape,
    /// so equal checksums after every phase mean equal layouts: leaf size classes, bursts
    /// into branch nodes, and merges back into leaves on erase. `find_common` returns the
    /// same entry (key and value) as C++.
    #[test]
    fn differential_vs_cpp_highs_hash_tree() {
        let expected = include_str!("../../tests/fixtures/mip/unit/hash_tree.out");
        let mut got = String::new();
        let nxt = |st: &mut u64| {
            *st ^= *st << 13;
            *st ^= *st >> 7;
            *st ^= *st << 17;
            *st
        };

        // Scenario A: random insert/erase with grow and shrink phases.
        let mut st: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut t: HighsHashTree<HighsInt, HighsInt> = HighsHashTree::default();
        for step in 0..60000i32 {
            let r = nxt(&mut st);
            let key = (nxt(&mut st) % 3000) as HighsInt;
            let phase = (step / 10000) % 2;
            let ins = if phase == 0 { r % 4 != 0 } else { r % 4 == 0 };
            if ins {
                t.insert(HighsHashTableEntry::new(key, step));
            } else {
                t.erase(&key);
            }
            if step % 2000 == 1999 {
                let (n, h) = fnv_of(&t);
                got += &format!("A {step} {n} {h}\n");
            }
        }

        // Scenario B: small key-only tree, full iteration order.
        st = 12345;
        let mut s: HighsHashTree<HighsInt, ()> = HighsHashTree::default();
        for _ in 0..120 {
            s.insert(HighsHashTableEntry::new(
                (nxt(&mut st) % 100000) as HighsInt,
                (),
            ));
        }
        for _ in 0..40 {
            s.erase(&((nxt(&mut st) % 100000) as HighsInt));
        }
        got += "B";
        s.for_each(&mut |k, _| got += &format!(" {k}"));
        got += "\n";

        // Scenario C: find_common between random trees.
        st = 777;
        let mut fc = Fnv::new();
        for trial in 0..400 {
            let mut t1: HighsHashTree<HighsInt, HighsInt> = HighsHashTree::default();
            let mut t2: HighsHashTree<HighsInt, HighsInt> = HighsHashTree::default();
            let n1 = 1 + nxt(&mut st) % 200;
            let n2 = 1 + nxt(&mut st) % 200;
            for i in 0..n1 {
                t1.insert(HighsHashTableEntry::new(
                    (nxt(&mut st) % 5000) as HighsInt,
                    i as i32,
                ));
            }
            for i in 0..n2 {
                t2.insert(HighsHashTableEntry::new(
                    (nxt(&mut st) % 5000) as HighsInt,
                    -(i as i32),
                ));
            }
            let out: i64 = match t1.find_common(&t2, false) {
                Some(e) => i64::from(*e.key()) * 1000 + i64::from(*e.value()),
                None => -1,
            };
            fc.add(out as u64);
            if trial < 10 {
                got += &format!("C {trial} {out}\n");
            }
        }
        got += &format!("Cfnv {}\n", fc.0);

        // Scenario D: fill, then drain in a pseudo-random order.
        st = 4242;
        let mut d: HighsHashTree<HighsInt, HighsInt> = HighsHashTree::default();
        let mut keys: Vec<HighsInt> = Vec::new();
        for i in 0..1500 {
            let k = (nxt(&mut st) % 1_000_000) as HighsInt;
            if d.insert(HighsHashTableEntry::new(k, i)) {
                keys.push(k);
            }
        }
        let mut i = keys.len();
        while i > 1 {
            let j = (nxt(&mut st) % i as u64) as usize;
            keys.swap(i - 1, j);
            i -= 1;
        }
        for i in 0..keys.len() {
            d.erase(&keys[i]);
            let left = keys.len() - i - 1;
            if left.is_multiple_of(50) || left < 60 {
                let (n, h) = fnv_of(&d);
                got += &format!("D {left} {n} {h}\n");
            }
        }
        got += &format!("Dempty {}\n", i32::from(d.empty()));

        for (lineno, (g, e)) in got.lines().zip(expected.lines()).enumerate() {
            assert_eq!(g, e, "line {}", lineno + 1);
        }
        assert_eq!(got.lines().count(), expected.lines().count());
    }

    #[test]
    fn for_each_mut_visits_in_for_each_order() {
        let mut tree: HighsHashTree<HighsInt, HighsInt> = HighsHashTree::default();
        for k in 0..500 {
            tree.insert(HighsHashTableEntry::new(k * 7, k));
        }
        let mut order = Vec::new();
        tree.for_each(&mut |k, _| order.push(*k));
        let mut order_mut = Vec::new();
        tree.for_each_mut(&mut |k, v| {
            order_mut.push(*k);
            *v += 1;
        });
        assert_eq!(order, order_mut);
        for k in 0..500 {
            assert_eq!(tree.find(&(k * 7)), Some(&(k + 1)));
        }
    }

    /// The C++ list-leaf erase (HighsHashTree.h:966-995) leaves a matching *last* node linked
    /// and only decrements `count`. Built by hand: list leaves need a 54-bit hash collision.
    #[test]
    fn list_leaf_erase_of_last_node_keeps_it_linked() {
        type T = HighsHashTree<HighsInt, HighsInt>;
        let mut node: HighsHashTreeNode<HighsInt, HighsInt> =
            HighsHashTreeNode::ListLeaf(Box::new(HighsHashTree_ListLeaf {
                first: HighsHashTree_ListNode {
                    next: Some(Box::new(HighsHashTree_ListNode {
                        next: None,
                        entry: HighsHashTableEntry::new(2, 20),
                    })),
                    entry: HighsHashTableEntry::new(1, 10),
                },
                count: 2,
            }));
        T::erase_recurse(&mut node, 0, kMaxDepth, &2);
        match &node {
            HighsHashTreeNode::ListLeaf(l) => assert_eq!(l.count, 1),
            _ => panic!("list leaf expected"),
        }
        assert_eq!(T::find_recurse(&node, 0, kMaxDepth, &2), Some(&20));
        // erasing the first node moves the (stale) second node into it; count hits 0
        T::erase_recurse(&mut node, 0, kMaxDepth, &1);
        assert_eq!(node.get_type(), kEmpty);
    }

    #[test]
    fn insert_or_get_returns_existing_value() {
        let mut tree: HighsHashTree<HighsInt, HighsInt> = HighsHashTree::default();
        let (v, inserted) = tree.insert_or_get(HighsHashTableEntry::new(7, 1));
        assert!(inserted);
        assert_eq!(*v, 1);
        let (v, inserted) = tree.insert_or_get(HighsHashTableEntry::new(7, 2));
        assert!(!inserted);
        assert_eq!(*v, 1);
    }

    /// HiGHS 1.16 #3327 (h16_hash_common): keys 57501, 91344 and 159347 share the top 16 bits of
    /// their hash (PR #3327 test `Highs_HashTree_findCommon_hash_collision`, diff lines 4-50).
    /// tree1 = {A, B}, tree2 = {B, C}; B is the only common key.
    fn h16_hash_common_trees() -> (
        HighsHashTree<HighsInt, HighsInt>,
        HighsHashTree<HighsInt, HighsInt>,
    ) {
        const KEY_A: HighsInt = 57501;
        const KEY_B: HighsInt = 91344;
        const KEY_C: HighsInt = 159347;
        let h = |k: HighsInt| HighsHashTree::<HighsInt, HighsInt>::compute_hash(&k) >> 48;
        assert_eq!(h(KEY_A), h(KEY_B));
        assert_eq!(h(KEY_A), h(KEY_C));
        let mut tree1: HighsHashTree<HighsInt, HighsInt> = HighsHashTree::default();
        let mut tree2: HighsHashTree<HighsInt, HighsInt> = HighsHashTree::default();
        tree1.insert(HighsHashTableEntry::new(KEY_A, 1));
        tree1.insert(HighsHashTableEntry::new(KEY_B, 2));
        tree2.insert(HighsHashTableEntry::new(KEY_B, 3));
        tree2.insert(HighsHashTableEntry::new(KEY_C, 4));
        (tree1, tree2)
    }

    #[test]
    fn h16_hash_common_on_finds_common_key_with_colliding_hashes() {
        let (tree1, tree2) = h16_hash_common_trees();
        assert!(tree1.contains(&91344) && tree2.contains(&91344));
        let found = tree1.find_common(&tree2, true).expect("common key 91344");
        assert_eq!(*found.key(), 91344);
        let found = tree2.find_common(&tree1, true).expect("common key 91344");
        assert_eq!(*found.key(), 91344);
    }

    #[test]
    fn h16_hash_common_off_matches_1_15_1_result() {
        // With the option off the 1.15.1 findCommonInLeaf runs unchanged on the same trees.
        let (tree1, tree2) = h16_hash_common_trees();
        assert!(tree1.find_common(&tree2, false).is_none());
        assert!(tree2.find_common(&tree1, false).is_none());
    }
}
