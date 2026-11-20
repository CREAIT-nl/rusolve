//! HiGHS source: util/HighsHash.h (HighsHashTableEntry, HighsHashTable)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Robin Hood open-addressing hash table with 8-bit metadata per slot (high
//! bit = occupied, low 7 bits = hash bits used to compute the distance from
//! the ideal slot). The slot layout decides iteration order, so the port
//! keeps the exact hash ([`HighsHashKey::highs_hash`] = C++
//! `HighsHashHelpers::hash`), the 7/8 growth rule, the 128-slot minimum,
//! backward-shift erase and the shrink-at-1/4 rule.
//!
//! Rust mapping:
//! * `std::unique_ptr<Entry> entries` (uninitialised storage) -> `entries:
//!   Vec<HighsHashTableEntry<K, V>>` of `capacity` default entries; a slot is
//!   live iff `occupied(metadata[pos])`. Placement-new / `~Entry()` become
//!   plain assignment (`entries[pos] = entry`), `swap(entry, entryArray[pos])`
//!   is `std::mem::swap(&mut entry, &mut self.entries[pos])`.
//! * `u64` arithmetic on positions uses `wrapping_sub` / `wrapping_add` where
//!   the C++ relies on unsigned wrap (`(pos - startPos) & tableSizeMask`,
//!   `(pos - metadata[pos]) & 0x7f`).
//! * `ValueType* find(key)` -> [`HighsHashTable::find`] (`Option<&V>`) and
//!   [`HighsHashTable::find_mut`] (`Option<&mut V>`); `nullptr` is `None`.
//! * `ValueType& operator[](key)` -> [`HighsHashTable::index_or_insert`]
//!   (`&mut V`). Internally it works with slot indices and borrows
//!   `self.entries[pos].value` only at the `return`.
//! * `insert(Args&&...)` -> [`HighsHashTable::insert`] taking an `Entry`
//!   (`HighsHashTableEntry::new(key, value)`), and a private
//!   `insert_entry` used by grow/shrink.
//! * `begin()/end()` iteration -> [`HighsHashTable::iter`] (slot order).
//!
//! LP presolve only instantiates `HighsHashTable<HighsInt, HighsInt>`
//! (`checkFillin`'s `fillinCache`, `detectParallelRowsAndCols`'
//! `numRowSingletons`) and never iterates it.

use crate::prelude::*;

/// Keys usable in a [`HighsHashTable`]: C++ `HighsHashHelpers::hash(key)`
/// and `HighsHashHelpers::equal` (memcmp, i.e. `PartialEq` for plain data).
pub trait HighsHashKey: Clone + Default + PartialEq {
    /// C++ `HighsHashHelpers::hash(const T&)` for this key type.
    fn highs_hash(&self) -> u64;
}

impl HighsHashKey for HighsInt {
    fn highs_hash(&self) -> u64 {
        hash_i32(*self)
    }
}

impl HighsHashKey for (u32, u32) {
    fn highs_hash(&self) -> u64 {
        hash_u32_pair(self.0, self.1)
    }
}

/// HiGHS: HighsHash.h::HighsHashTableEntry<K, V>
#[derive(Clone, Copy, Debug, Default)]
pub struct HighsHashTableEntry<K, V> {
    pub key_: K,
    pub value_: V,
}

impl<K, V> HighsHashTableEntry<K, V> {
    pub fn new(key: K, value: V) -> Self {
        Self {
            key_: key,
            value_: value,
        }
    }
    /// HiGHS: HighsHashTableEntry::key
    pub fn key(&self) -> &K {
        &self.key_
    }
    /// HiGHS: HighsHashTableEntry::value
    pub fn value(&self) -> &V {
        &self.value_
    }
}

/// HiGHS: HighsHash.h::HighsHashTable<K, V>
///
/// C++ initialisers: `numElements = 0`; the default constructor calls
/// `makeEmptyTable(128)` (so use [`HighsHashTable::new`], which `Default`
/// forwards to). `tableSizeMask`/`numHashShift` are always set by
/// `makeEmptyTable`.
#[derive(Clone, Debug)]
pub struct HighsHashTable<K, V> {
    entries: Vec<HighsHashTableEntry<K, V>>,
    metadata: Vec<u8>,
    table_size_mask: u64,
    num_hash_shift: u64,
    num_elements: u64,
}

impl<K: HighsHashKey, V: Clone + Default> Default for HighsHashTable<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: HighsHashKey, V: Clone + Default> HighsHashTable<K, V> {
    /// HiGHS: HighsHash.h::HighsHashTable::HighsHashTable() (line 986): `makeEmptyTable(128)`.
    pub fn new() -> Self {
        let mut t = Self {
            entries: Vec::new(),
            metadata: Vec::new(),
            table_size_mask: 0,
            num_hash_shift: 0,
            num_elements: 0,
        };
        t.make_empty_table(128);
        t
    }

    /// HiGHS: HighsHash.h::HighsHashTable::HighsHashTable(u64 minCapacity) (lines 987-991, 5 lines)
    ///
    /// `std::ceil(std::log2(std::max(128.0, 8 * double(minCapacity) / 7)))`:
    /// use `f64::log2`/`f64::ceil` and `cpp_max`; keep the evaluation order
    /// `(8.0 * min_capacity as f64) / 7.0`.
    pub fn with_min_capacity(min_capacity: u64) -> Self {
        let init_capacity: u64 = 1u64
            << (cpp_max(128.0, (8.0 * min_capacity as f64) / 7.0)
                .log2()
                .ceil() as u64);
        let mut t = Self {
            entries: Vec::new(),
            metadata: Vec::new(),
            table_size_mask: 0,
            num_hash_shift: 0,
            num_elements: 0,
        };
        t.make_empty_table(init_capacity);
        t
    }

    /// HiGHS: HighsHash.h::HighsHashTable::toMetadata (line 1041)
    fn to_metadata(&self, hash: u64) -> u8 {
        ((hash >> self.num_hash_shift) | 0x80) as u8
    }

    /// HiGHS: HighsHash.h::HighsHashTable::maxDistance (line 1043)
    const fn max_distance() -> u64 {
        127
    }

    /// HiGHS: HighsHash.h::HighsHashTable::makeEmptyTable (lines 1045-1054, 10 lines)
    ///
    /// `numHashShift = 64 - HighsHashHelpers::log2i(capacity)` uses
    /// `log2i_u64`. `new u8[capacity]{}` -> `vec![0u8; capacity]`; entries
    /// -> `vec![HighsHashTableEntry::default(); capacity]`.
    fn make_empty_table(&mut self, capacity: u64) {
        self.table_size_mask = capacity - 1;
        self.num_hash_shift = (64 - log2i_u64(capacity)) as u64;
        debug_assert!(capacity == (1u64 << (64 - self.num_hash_shift)));
        self.num_elements = 0;

        self.metadata = vec![0u8; capacity as usize];
        self.entries = vec![HighsHashTableEntry::default(); capacity as usize];
    }

    /// HiGHS: HighsHash.h::HighsHashTable::occupied (line 1056)
    fn occupied(meta: u8) -> bool {
        meta & 0x80 != 0
    }

    /// HiGHS: HighsHash.h::HighsHashTable::distanceFromIdealSlot (lines 1058-1075, 18 lines)
    fn distance_from_ideal_slot(&self, pos: u64) -> u64 {
        pos.wrapping_sub(self.metadata[pos as usize] as u64) & 0x7f
    }

    /// HiGHS: HighsHash.h::HighsHashTable::growTable (lines 1077-1086, 10 lines)
    ///
    /// `insert(std::move(oldEntries[i]))` -> `self.insert_entry(entry)`.
    fn grow_table(&mut self) {
        let mut old_entries = std::mem::take(&mut self.entries);
        let old_metadata = std::mem::take(&mut self.metadata);
        let old_capactiy = self.table_size_mask + 1;

        self.make_empty_table(2 * old_capactiy);

        for i in 0..old_capactiy as usize {
            if Self::occupied(old_metadata[i]) {
                let entry = std::mem::take(&mut old_entries[i]);
                self.insert_entry(entry);
            }
        }
    }

    /// HiGHS: HighsHash.h::HighsHashTable::shrinkTable (lines 1088-1097, 10 lines)
    fn shrink_table(&mut self) {
        let mut old_entries = std::mem::take(&mut self.entries);
        let old_metadata = std::mem::take(&mut self.metadata);
        let old_capactiy = self.table_size_mask + 1;

        self.make_empty_table(old_capactiy / 2);

        for i in 0..old_capactiy as usize {
            if Self::occupied(old_metadata[i]) {
                let entry = std::mem::take(&mut old_entries[i]);
                self.insert_entry(entry);
            }
        }
    }

    /// HiGHS: HighsHash.h::HighsHashTable::findPosition (lines 1099-1121, 23 lines)
    ///
    /// Out-params `meta`, `startPos`, `maxPos`, `pos` are returned as the
    /// tuple `(found, meta, start_pos, max_pos, pos)`.
    #[inline]
    fn find_position(&self, key: &K) -> (bool, u8, u64, u64, u64) {
        let hash = key.highs_hash();
        let start_pos = hash >> self.num_hash_shift;
        let max_pos = (start_pos + Self::max_distance()) & self.table_size_mask;
        let meta = self.to_metadata(hash);

        // PORT NOTE (perf): the slot arrays and the mask are read once into locals; each
        // slot's metadata byte is loaded once and `distanceFromIdealSlot` is computed from it.
        let mask = self.table_size_mask;
        let metadata = &self.metadata[..];
        let entries = &self.entries[..];
        let mut pos = start_pos;
        loop {
            let m = metadata[pos as usize];
            if !Self::occupied(m) {
                return (false, meta, start_pos, max_pos, pos);
            }
            if m == meta && key == entries[pos as usize].key() {
                return (true, meta, start_pos, max_pos, pos);
            }

            let current_distance = pos.wrapping_sub(start_pos) & mask;
            if current_distance > (pos.wrapping_sub(m as u64) & 0x7f) {
                return (false, meta, start_pos, max_pos, pos);
            }

            pos = (pos + 1) & mask;
            if pos == max_pos {
                break;
            }
        }
        (false, meta, start_pos, max_pos, pos)
    }

    /// HiGHS: HighsHash.h::HighsHashTable::clear (lines 1124-1137, 14 lines)
    ///
    /// Destructor calls are no-ops; `std::memset(metadata, 0, 128)` ->
    /// `self.metadata.fill(0)`.
    pub fn clear(&mut self) {
        if self.num_elements != 0 {
            let capacity = self.table_size_mask + 1;
            if capacity == 128 {
                self.metadata.fill(0);
                self.num_elements = 0;
            } else {
                self.make_empty_table(128);
            }
        }
    }

    /// HiGHS: HighsHash.h::HighsHashTable::find (const) (lines 1139-1146, 8 lines)
    #[inline]
    pub fn find(&self, key: &K) -> Option<&V> {
        let (found, _meta, _start_pos, _max_pos, pos) = self.find_position(key);
        if found {
            Some(&self.entries[pos as usize].value_)
        } else {
            None
        }
    }

    /// HiGHS: HighsHash.h::HighsHashTable::find (non-const) (lines 1148-1155, 8 lines)
    pub fn find_mut(&mut self, key: &K) -> Option<&mut V> {
        let (found, _meta, _start_pos, _max_pos, pos) = self.find_position(key);
        if found {
            Some(&mut self.entries[pos as usize].value_)
        } else {
            None
        }
    }

    /// HiGHS: HighsHash.h::HighsHashTable::operator[] (lines 1157-1197, 41 lines)
    ///
    /// `ValueType& insertLocation = entryArray[pos].value()` is a slot index
    /// `insert_location = pos`, returned as `&mut self.entries[pos].value_`.
    /// `Entry entry(key)` value-initialises the value: `V::default()`.
    /// Recursive `return (*this)[key];` -> `return self.index_or_insert(key);`.
    pub fn index_or_insert(&mut self, key: &K) -> &mut V {
        let (found, mut meta, mut start_pos, mut max_pos, mut pos) = self.find_position(key);
        if found {
            return &mut self.entries[pos as usize].value_;
        }

        if self.num_elements == ((self.table_size_mask + 1) * 7) / 8 || pos == max_pos {
            self.grow_table();
            return self.index_or_insert(key);
        }

        let insert_location = pos;
        let mut entry = HighsHashTableEntry::new(key.clone(), V::default());
        self.num_elements += 1;

        loop {
            if !Self::occupied(self.metadata[pos as usize]) {
                self.metadata[pos as usize] = meta;
                self.entries[pos as usize] = entry;
                return &mut self.entries[insert_location as usize].value_;
            }

            let current_distance = pos.wrapping_sub(start_pos) & self.table_size_mask;
            let distance_of_current_occupant = self.distance_from_ideal_slot(pos);
            if current_distance > distance_of_current_occupant {
                // steal the position
                std::mem::swap(&mut entry, &mut self.entries[pos as usize]);
                std::mem::swap(&mut meta, &mut self.metadata[pos as usize]);

                start_pos = pos.wrapping_sub(distance_of_current_occupant) & self.table_size_mask;
                max_pos = (start_pos + Self::max_distance()) & self.table_size_mask;
            }
            pos = (pos + 1) & self.table_size_mask;
            if pos == max_pos {
                break;
            }
        }

        self.grow_table();
        self.insert(entry);
        self.index_or_insert(key)
    }

    /// HiGHS: HighsHash.h::HighsHashTable::insert (lines 1199-1239, 41 lines)
    ///
    /// The variadic `Entry entry(std::forward<Args>(args)...)` is the
    /// caller-built `entry`. Recursive `insert(std::move(entry))` calls stay
    /// recursive.
    pub fn insert(&mut self, mut entry: HighsHashTableEntry<K, V>) -> bool {
        let (found, mut meta, mut start_pos, mut max_pos, mut pos) =
            self.find_position(entry.key());
        if found {
            return false;
        }

        if self.num_elements == ((self.table_size_mask + 1) * 7) / 8 || pos == max_pos {
            self.grow_table();
            return self.insert(entry);
        }

        self.num_elements += 1;

        loop {
            if !Self::occupied(self.metadata[pos as usize]) {
                self.metadata[pos as usize] = meta;
                self.entries[pos as usize] = entry;
                return true;
            }

            let current_distance = pos.wrapping_sub(start_pos) & self.table_size_mask;
            let distance_of_current_occupant = self.distance_from_ideal_slot(pos);
            if current_distance > distance_of_current_occupant {
                // steal the position
                std::mem::swap(&mut entry, &mut self.entries[pos as usize]);
                std::mem::swap(&mut meta, &mut self.metadata[pos as usize]);

                start_pos = pos.wrapping_sub(distance_of_current_occupant) & self.table_size_mask;
                max_pos = (start_pos + Self::max_distance()) & self.table_size_mask;
            }
            pos = (pos + 1) & self.table_size_mask;
            if pos == max_pos {
                break;
            }
        }

        self.grow_table();
        self.insert(entry);
        true
    }

    /// Internal alias used by `grow_table` / `shrink_table` (C++ calls the
    /// same `insert`); forwards to [`Self::insert`].
    fn insert_entry(&mut self, entry: HighsHashTableEntry<K, V>) -> bool {
        self.insert(entry)
    }

    /// HiGHS: HighsHash.h::HighsHashTable::erase (lines 1241-1272, 32 lines)
    ///
    /// `entryArray[pos] = std::move(entryArray[shift])` ->
    /// `self.entries[pos] = std::mem::take(&mut self.entries[shift])`.
    pub fn erase(&mut self, key: &K) -> bool {
        let (found, _meta, _start_pos, _max_pos, mut pos) = self.find_position(key);
        if !found {
            return false;
        }
        // delete element at position pos: `~Entry()` is a no-op here, the
        // slot's stale value is never read because its metadata is cleared.
        self.metadata[pos as usize] = 0;

        // retain at least a quarter of slots occupied, otherwise shrink the table
        // if its not at its minimum size already
        self.num_elements -= 1;
        let capacity = self.table_size_mask + 1;
        if capacity != 128 && self.num_elements < capacity / 4 {
            self.shrink_table();
            return true;
        }

        // shift elements after pos backwards
        loop {
            let shift = (pos + 1) & self.table_size_mask;
            if !Self::occupied(self.metadata[shift as usize]) {
                return true;
            }

            let dist = self.distance_from_ideal_slot(shift);
            if dist == 0 {
                return true;
            }

            self.entries[pos as usize] = std::mem::take(&mut self.entries[shift as usize]);
            self.metadata[pos as usize] = self.metadata[shift as usize];
            self.metadata[shift as usize] = 0;
            pos = shift;
        }
    }

    /// HiGHS: HighsHash.h::HighsHashTable::size (line 1274)
    pub fn size(&self) -> i64 {
        self.num_elements as i64
    }

    /// HiGHS: HighsHash.h::HighsHashTable::begin/end (lines 1008-1038):
    /// occupied slots in slot order.
    pub fn iter(&self) -> impl Iterator<Item = &HighsHashTableEntry<K, V>> + '_ {
        let n = if self.num_elements == 0 {
            0
        } else {
            self.metadata.len()
        };
        (0..n)
            .filter(move |&i| Self::occupied(self.metadata[i]))
            .map(move |i| &self.entries[i])
    }
}
