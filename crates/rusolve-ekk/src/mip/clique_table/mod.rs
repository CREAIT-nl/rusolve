//! HiGHS source: mip/HighsCliqueTable.h (class definitions and header-inline methods)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

mod part1;
pub use part1::*;
mod part2;
pub use part2::*;

/// HiGHS: HighsCliqueTable.h::HighsCliqueTable (lines 29-338)
#[derive(Clone, Default, Debug)]
pub struct HighsCliqueTable {
    /// C++: `std::vector<CliqueVar> cliqueentries`
    pub cliqueentries: Vec<HighsCliqueTable_CliqueVar>,
    /// C++: `std::vector<HighsHashTree<HighsInt, HighsInt>> invertedHashList`
    pub inverted_hash_list: Vec<HighsHashTree<HighsInt, HighsInt>>,
    /// C++: `std::vector<HighsHashTree<HighsInt>> invertedHashListSizeTwo`
    pub inverted_hash_list_size_two: Vec<HighsHashTree<HighsInt, ()>>,
    /// C++: `HighsHashTable<std::pair<CliqueVar, CliqueVar>, HighsInt> sizeTwoCliques`
    /// PORT NOTE (perf): keyed by the raw bitfield words of the two vars (see `sorted_edge`).
    pub size_two_cliques: HighsHashTable<(u32, u32), HighsInt>,
    /// HiGHS 1.16 (#3540): HighsCliqueTable.h::HighsCliqueTable (lines 79-86 of the 1.16.0 header)
    /// C++: `std::vector<SizeTwoLinks> sizeTwoLinks` (inert unless `h16_clique_links`)
    pub size_two_links: Vec<HighsCliqueTable_SizeTwoLinks>,
    /// C++: `std::set<std::pair<HighsInt, HighsInt>> freespaces`
    pub freespaces: BTreeSet<(HighsInt, HighsInt)>,
    /// C++: `std::vector<HighsInt> freeslots`
    pub freeslots: Vec<HighsInt>,
    /// C++: `std::vector<Clique> cliques`
    pub cliques: Vec<HighsCliqueTable_Clique>,
    /// C++: `std::vector<HighsInt> numcliquesvar`
    pub numcliquesvar: Vec<HighsInt>,
    /// C++: `std::vector<CliqueVar> infeasvertexstack`
    pub infeasvertexstack: Vec<HighsCliqueTable_CliqueVar>,
    /// C++: `std::vector<HighsInt> colsubstituted`
    pub colsubstituted: Vec<HighsInt>,
    /// C++: `std::vector<Substitution> substitutions`
    pub substitutions: Vec<HighsCliqueTable_Substitution>,
    /// C++: `std::vector<HighsInt> deletedrows`
    pub deletedrows: Vec<HighsInt>,
    /// C++: `std::vector<std::pair<HighsInt, CliqueVar>> cliqueextensions`
    pub cliqueextensions: Vec<(HighsInt, HighsCliqueTable_CliqueVar)>,
    /// C++: `std::vector<uint8_t> iscandidate`
    pub iscandidate: Vec<u8>,
    /// C++: `std::vector<uint8_t> colDeleted`
    pub col_deleted: Vec<u8>,
    /// C++: `std::vector<uint32_t> cliquehits`
    pub cliquehits: Vec<u32>,
    /// C++: `std::vector<HighsInt> cliquehitinds`
    pub cliquehitinds: Vec<HighsInt>,
    /// C++: `HighsRandom randgen`
    pub randgen: HighsRandom,
    /// C++: `HighsInt nfixings`
    pub nfixings: HighsInt,
    /// C++: `HighsInt numEntries`
    pub num_entries: HighsInt,
    /// C++: `HighsInt maxEntries`
    pub max_entries: HighsInt,
    /// C++: `HighsInt minEntriesForParallelism`
    pub min_entries_for_parallelism: HighsInt,
    /// C++: `bool inPresolve`
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.h `bool inPresolveProbing` (line 102):
    /// with `mip_presolve_clique_table` it is set only around presolve probing and
    /// solution enumeration (CLQP-E); the field keeps its 1.15.1 name.
    pub in_presolve: bool,
    /// C++: `bool allowParallel`
    pub allow_parallel: bool,
    /// C++: `int64_t numNeighbourhoodQueries`
    pub num_neighbourhood_queries: i64,
    /// HiGHS 1.16 #3540 (h16_clique_links): fixed at construction by `set_h16_flags` and copied by
    /// `copy_h16_flags_from`; false reproduces 1.15.1
    pub h16_clique_links: bool,
    /// HiGHS 1.16 #3327 (h16_hash_common): fixed at construction by `set_h16_flags` and copied by
    /// `copy_h16_flags_from`; false reproduces 1.15.1 (`find_common` without the 1.16.0 leaf compare)
    pub h16_hash_common: bool,
    /// HiGHS #3312 (mip_presolve_clique_table): fixed at construction by `set_h16_flags` and
    /// copied by `copy_h16_flags_from`; false reproduces 1.15.1. When true, `in_presolve` has the
    /// 1.16.0 meaning `inPresolveProbing` (set only around presolve probing and enumeration).
    pub x1b_presolve_clique_table: bool,
}

/// HiGHS 1.16 (#3540): HighsCliqueTable.h::HighsCliqueTable::SizeTwoLinks (lines 82-85 of the 1.16.0 header)
/// C++: `struct SizeTwoLinks { HighsInt prev; HighsInt next; }`
#[derive(Clone, Copy, Default, Debug)]
pub struct HighsCliqueTable_SizeTwoLinks {
    /// C++: `HighsInt prev`
    pub prev: HighsInt,
    /// C++: `HighsInt next`
    pub next: HighsInt,
}

/// HiGHS: HighsCliqueTable.h::HighsCliqueTable::CliqueVar (lines 31-54)
///
/// C++ is the 4-byte bitfield `HighsUInt col : 31; HighsUInt val : 1;`. The port stores
/// the same object representation in one `u32` (little-endian GCC/Clang layout): `col`
/// in the low 31 bits, `val` in the top bit. Construction truncates exactly as the C++
/// bitfield assignment does (`col & 0x7FFF_FFFF`, `val & 1`).
#[derive(Clone, Copy, Default)]
pub struct HighsCliqueTable_CliqueVar {
    /// C++ bitfields `col : 31` (bits 0-30) and `val : 1` (bit 31).
    bits: u32,
}

impl std::fmt::Debug for HighsCliqueTable_CliqueVar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HighsCliqueTable_CliqueVar")
            .field("col", &self.col())
            .field("val", &self.val())
            .finish()
    }
}

impl PartialEq for HighsCliqueTable_CliqueVar {
    fn eq(&self, other: &Self) -> bool {
        self.op_eq(other)
    }
}

/// HiGHS: HighsCliqueTable.h::HighsCliqueTable::Clique (lines 55-62)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct HighsCliqueTable_Clique {
    /// C++: `HighsInt start`
    pub start: HighsInt,
    /// C++: `HighsInt end`
    pub end: HighsInt,
    /// C++: `HighsInt origin`
    pub origin: HighsInt,
    /// C++: `HighsInt numZeroFixed`
    pub num_zero_fixed: HighsInt,
    /// C++: `bool equality`
    pub equality: bool,
}

/// HiGHS: HighsCliqueTable.h::HighsCliqueTable::Substitution (lines 64-67)
#[derive(Clone, Default, Debug)]
pub struct HighsCliqueTable_Substitution {
    /// C++: `HighsInt substcol`
    pub substcol: HighsInt,
    /// C++: `CliqueVar replace`
    pub replace: HighsCliqueTable_CliqueVar,
}

/// HiGHS: HighsCliqueTable.h::HighsCliqueTable::BronKerboschData (lines 121-144)
#[derive(Clone, Debug)]
pub struct HighsCliqueTable_BronKerboschData<'a> {
    /// C++ `const std::vector<double>& sol`.
    pub sol: &'a [f64],
    /// C++: `std::vector<CliqueVar> P`
    pub p: Vec<HighsCliqueTable_CliqueVar>,
    /// C++: `std::vector<CliqueVar> R`
    pub r: Vec<HighsCliqueTable_CliqueVar>,
    /// C++: `std::vector<CliqueVar> Z`
    pub z: Vec<HighsCliqueTable_CliqueVar>,
    /// C++: `std::vector<std::vector<CliqueVar>> cliques`
    pub cliques: Vec<Vec<HighsCliqueTable_CliqueVar>>,
    /// C++: `std::vector<HighsInt> neighbourhoodInds`
    pub neighbourhood_inds: Vec<HighsInt>,
    /// C++: `double wR = 0.0`
    pub w_r: f64,
    /// C++: `double minW = 1.05`
    pub min_w: f64,
    /// C++: `double feastol = 1e-6`
    pub feastol: f64,
    /// C++: `HighsInt ncalls = 0`
    pub ncalls: HighsInt,
    /// C++: `HighsInt maxcalls = 10000`
    pub maxcalls: HighsInt,
    /// C++: `HighsInt maxcliques = 100`
    pub maxcliques: HighsInt,
    /// C++: `int64_t maxNeighbourhoodQueries = std :: numeric_limits < int64_t > :: max ( )`
    pub max_neighbourhood_queries: i64,
    /// C++: `int64_t numNeighbourhoodQueries = 0`
    pub num_neighbourhood_queries: i64,
}

impl Default for HighsCliqueTable_BronKerboschData<'_> {
    fn default() -> Self {
        Self {
            sol: Default::default(),
            p: Default::default(),
            r: Default::default(),
            z: Default::default(),
            cliques: Default::default(),
            neighbourhood_inds: Default::default(),
            w_r: 0.0,
            min_w: 1.05,
            feastol: 1e-6,
            ncalls: 0,
            maxcalls: 10000,
            maxcliques: 100,
            max_neighbourhood_queries: i64::MAX,
            num_neighbourhood_queries: 0,
        }
    }
}

/// HiGHS: HighsCliqueTable.cpp::ThreadNeighbourhoodQueryData (lines 475-478)
#[derive(Clone, Default, Debug)]
pub struct ThreadNeighbourhoodQueryData {
    /// C++: `int64_t numQueries`
    pub num_queries: i64,
    /// C++: `std::vector<HighsInt> neighbourhoodInds`
    pub neighbourhood_inds: Vec<HighsInt>,
}

impl HighsCliqueTable_CliqueVar {
    /// Bitfield constructor from already-unsigned fields: `col : 31` and `val : 1` keep
    /// the low 31 bits and the low bit, as the C++ member initialisers do.
    #[inline(always)]
    pub const fn from_parts(col: HighsUInt, val: HighsUInt) -> Self {
        Self {
            bits: (col & 0x7FFF_FFFF) | ((val & 1) << 31),
        }
    }

    /// C++ bitfield read `col` (31 bits, zero-extended).
    #[inline(always)]
    pub const fn col(&self) -> HighsUInt {
        self.bits & 0x7FFF_FFFF
    }

    /// C++ bitfield read `val` (1 bit).
    #[inline(always)]
    pub const fn val(&self) -> HighsUInt {
        self.bits >> 31
    }

    /// C++ bitfield write `col = c` (truncated to 31 bits; `val` unchanged).
    #[inline(always)]
    pub fn set_col(&mut self, col: HighsUInt) {
        self.bits = (self.bits & 0x8000_0000) | (col & 0x7FFF_FFFF);
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::CliqueVar::index (lines 40-40)
    /// C++: `HighsInt index() const`
    #[inline(always)]
    pub fn index(&self) -> HighsInt {
        2 * (self.col() as HighsInt) + (self.val() as HighsInt)
    }

    /// C++ raw 32-bit word of the `CliqueVar` bitfield struct, where `col : 31` occupies
    /// the low bits and `val : 1` the top bit (little-endian GCC/Clang layout). This is the
    /// object representation `HighsHashHelpers::hash` and `memcmp` see.
    #[inline(always)]
    pub fn word(&self) -> u32 {
        self.bits
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::CliqueVar::weight (lines 42-44)
    /// C++: `double weight(const std::vector<double> & sol) const`
    #[inline]
    pub fn weight(&self, sol: &[f64]) -> f64 {
        if self.val() != 0 {
            sol[self.col() as usize]
        } else {
            1.0 - sol[self.col() as usize]
        }
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::CliqueVar::complement (lines 46-46)
    /// C++: `CliqueVar complement() const` = `CliqueVar(col, 1 - val)`; with `val` one bit
    /// wide this flips bit 31.
    #[inline(always)]
    pub fn complement(&self) -> HighsCliqueTable_CliqueVar {
        Self {
            bits: self.bits ^ 0x8000_0000,
        }
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::CliqueVar::operator== (lines 48-50)
    /// C++: `bool operator==(const CliqueVar & other) const`
    #[inline(always)]
    pub fn op_eq(&self, other: &HighsCliqueTable_CliqueVar) -> bool {
        self.index() == other.index()
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::CliqueVar::CliqueVar (lines 52-52)
    /// C++: `CliqueVar(HighsInt col, HighsInt val)`
    #[inline(always)]
    pub fn new(col: HighsInt, val: HighsInt) -> Self {
        Self::from_parts(col as HighsUInt, val as HighsUInt)
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::CliqueVar::CliqueVar (lines 53-53)
    /// C++: `CliqueVar()`
    pub fn new_2() -> Self {
        Self::default()
    }
}

impl HighsCliqueTable_Clique {
    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::Clique::numActive (lines 61-61)
    /// C++: `HighsInt numActive() const`
    pub fn num_active(&self) -> HighsInt {
        self.end - self.start - self.num_zero_fixed
    }
}

impl HighsCliqueTable {
    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::findCommonCliqueId (lines 108-110)
    /// C++: `HighsInt findCommonCliqueId(CliqueVar v1, CliqueVar v2)`
    pub fn find_common_clique_id_2(
        &mut self,
        v1: HighsCliqueTable_CliqueVar,
        v2: HighsCliqueTable_CliqueVar,
    ) -> HighsInt {
        // PORT NOTE: C++ findCommonCliqueId(v1, v2) (non-const) passes the member
        // counter numNeighbourhoodQueries by reference; copy in and out around the
        // const overload.
        let mut num_queries = self.num_neighbourhood_queries;
        let id = self.find_common_clique_id(&mut num_queries, v1, v2);
        self.num_neighbourhood_queries = num_queries;
        id
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::HighsCliqueTable (lines 166-179)
    /// C++: `HighsCliqueTable(HighsInt ncols)`
    pub fn new(ncols: HighsInt) -> Self {
        let n = ncols as usize;
        Self {
            // C++ invertedHashList.resize(2 * ncols), invertedHashListSizeTwo likewise
            inverted_hash_list: (0..2 * n).map(|_| HighsHashTree::default()).collect(),
            inverted_hash_list_size_two: (0..2 * n).map(|_| HighsHashTree::default()).collect(),
            numcliquesvar: vec![0; 2 * n],
            colsubstituted: vec![0; n],
            col_deleted: vec![0; n],
            nfixings: 0,
            num_neighbourhood_queries: 0,
            num_entries: 0,
            max_entries: kHighsIInf,
            min_entries_for_parallelism: kHighsIInf,
            in_presolve: false,
            allow_parallel: true,
            // PORT NOTE: C++ `HighsRandom randgen` is default-constructed with seed 0,
            // which is HighsRandom(0) (initialise(0)), not the derived Default state.
            randgen: HighsRandom::new(0),
            ..Default::default()
        }
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::setPresolveFlag (lines 181-181)
    /// C++: `void setPresolveFlag(bool inPresolve)`
    pub fn set_presolve_flag(&mut self, in_presolve: bool) {
        self.in_presolve = in_presolve;
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getPresolveFlag (lines 183-183)
    /// C++: `bool getPresolveFlag() const`
    pub fn get_presolve_flag(&self) -> bool {
        self.in_presolve
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getNumEntries (lines 185-185)
    /// C++: `HighsInt getNumEntries() const`
    pub fn get_num_entries(&self) -> HighsInt {
        self.num_entries
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getRandgen (lines 187-187)
    /// C++: `HighsRandom & getRandgen()`
    pub fn get_randgen(&mut self) -> &mut HighsRandom {
        &mut self.randgen
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getNumNeighbourhoodQueries (lines 189-189)
    /// C++: `int64_t & getNumNeighbourhoodQueries()`
    pub fn get_num_neighbourhood_queries(&mut self) -> &mut i64 {
        &mut self.num_neighbourhood_queries
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getDeletedRows (lines 217-217)
    /// C++: `std::vector<HighsInt> & getDeletedRows()`
    pub fn get_deleted_rows(&mut self) -> &mut Vec<HighsInt> {
        &mut self.deletedrows
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getDeletedRows (lines 219-219)
    /// C++: `const std::vector<HighsInt> & getDeletedRows() const`
    pub fn get_deleted_rows_2(&self) -> &[HighsInt] {
        &self.deletedrows
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getSubstitutions (lines 221-221)
    /// C++: `std::vector<Substitution> & getSubstitutions()`
    pub fn get_substitutions(&mut self) -> &mut Vec<HighsCliqueTable_Substitution> {
        &mut self.substitutions
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getSubstitutions (lines 223-225)
    /// C++: `const std::vector<Substitution> & getSubstitutions() const`
    pub fn get_substitutions_2(&self) -> &[HighsCliqueTable_Substitution] {
        &self.substitutions
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getSubstitution (lines 227-230)
    /// C++: `const Substitution * getSubstitution(HighsInt col) const`
    pub fn get_substitution(&self, col: HighsInt) -> Option<&HighsCliqueTable_Substitution> {
        let idx = self.colsubstituted[col as usize];
        if idx != 0 {
            Some(&self.substitutions[(idx - 1) as usize])
        } else {
            None
        }
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getCliqueExtensions (lines 232-234)
    /// C++: `std::vector<std::pair<HighsInt, CliqueVar>> & getCliqueExtensions()`
    pub fn get_clique_extensions(&mut self) -> &mut Vec<(HighsInt, HighsCliqueTable_CliqueVar)> {
        &mut self.cliqueextensions
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getCliqueExtensions (lines 236-239)
    /// C++: `const std::vector<std::pair<HighsInt, CliqueVar>> & getCliqueExtensions() const`
    pub fn get_clique_extensions_2(&self) -> &[(HighsInt, HighsCliqueTable_CliqueVar)] {
        &self.cliqueextensions
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::setMaxEntries (lines 241-243)
    /// C++: `void setMaxEntries(HighsInt numNz)`
    pub fn set_max_entries(&mut self, num_nz: HighsInt) {
        self.max_entries = 2000000 + 10 * num_nz;
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::setMinEntriesForParallelism (lines 245-247)
    /// C++: `void setMinEntriesForParallelism(HighsInt minEntriesForParallelism)`
    pub fn set_min_entries_for_parallelism(&mut self, min_entries_for_parallelism: HighsInt) {
        self.min_entries_for_parallelism = min_entries_for_parallelism;
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::isFull (lines 249-249)
    /// C++: `bool isFull() const`
    pub fn is_full(&self) -> bool {
        self.num_entries >= self.max_entries
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::getNumFixings (lines 251-251)
    /// C++: `HighsInt getNumFixings() const`
    pub fn get_num_fixings(&self) -> HighsInt {
        self.nfixings
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::haveCommonClique (lines 272-275)
    /// C++: `bool haveCommonClique(CliqueVar v1, CliqueVar v2)`
    pub fn have_common_clique(
        &mut self,
        v1: HighsCliqueTable_CliqueVar,
        v2: HighsCliqueTable_CliqueVar,
    ) -> bool {
        if v1.col() == v2.col() {
            return false;
        }
        self.find_common_clique_id_2(v1, v2) != -1
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::haveCommonClique (lines 277-280)
    /// C++: `bool haveCommonClique(int64_t & numQueries, CliqueVar v1, CliqueVar v2) const`
    pub fn have_common_clique_2(
        &self,
        num_queries: &mut i64,
        v1: HighsCliqueTable_CliqueVar,
        v2: HighsCliqueTable_CliqueVar,
    ) -> bool {
        if v1.col() == v2.col() {
            return false;
        }
        self.find_common_clique_id(num_queries, v1, v2) != -1
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::findCommonClique (lines 282-292)
    /// C++: `std::pair<const CliqueVar *, HighsInt> findCommonClique(CliqueVar v1, CliqueVar v2)`
    pub fn find_common_clique(
        &mut self,
        v1: HighsCliqueTable_CliqueVar,
        v2: HighsCliqueTable_CliqueVar,
    ) -> (&[HighsCliqueTable_CliqueVar], HighsInt) {
        if v1 == v2 {
            return (&[], 0);
        }
        let clq = self.find_common_clique_id_2(v1, v2);
        if clq == -1 {
            return (&[], 0);
        }
        let start = self.cliques[clq as usize].start as usize;
        let end = self.cliques[clq as usize].end as usize;
        (&self.cliqueentries[start..end], (end - start) as HighsInt)
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::numCliques (lines 327-327)
    /// C++: `HighsInt numCliques() const`
    pub fn num_cliques(&self) -> HighsInt {
        (self.cliques.len() - self.freeslots.len()) as HighsInt
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::numCliques (lines 329-329)
    /// C++: `HighsInt numCliques(CliqueVar v) const`
    pub fn num_cliques_2(&self, v: HighsCliqueTable_CliqueVar) -> HighsInt {
        self.numcliquesvar[v.index() as usize]
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::numCliques (lines 331-333)
    /// C++: `HighsInt numCliques(HighsInt col, bool val) const`
    pub fn num_cliques_3(&self, col: HighsInt, val: bool) -> HighsInt {
        self.numcliquesvar[HighsCliqueTable_CliqueVar::new(col, val as HighsInt).index() as usize]
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::setAllowParallel (lines 335-337)
    /// C++: `void setAllowParallel(const bool allowParallel)`
    pub fn set_allow_parallel(&mut self, allow_parallel: bool) {
        self.allow_parallel = allow_parallel;
    }
}

impl<'a> HighsCliqueTable_BronKerboschData<'a> {
    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::BronKerboschData::stop (lines 138-141)
    /// C++: `bool stop() const`
    pub fn stop(&self) -> bool {
        self.maxcalls == self.ncalls
            || self.cliques.len() as HighsInt == self.maxcliques
            || self.num_neighbourhood_queries > self.max_neighbourhood_queries
    }

    /// HiGHS: HighsCliqueTable.h::HighsCliqueTable::BronKerboschData::BronKerboschData (lines 143-143)
    /// C++: `BronKerboschData(const std::vector<double> & sol)`
    pub fn new(sol: &'a [f64]) -> Self {
        Self {
            sol,
            ..Self::default()
        }
    }
}

#[cfg(test)]
mod clique_var_layout_tests {
    use super::HighsCliqueTable_CliqueVar as CliqueVar;

    #[test]
    fn packed_like_the_cpp_bitfield() {
        // C++ `HighsUInt col : 31; HighsUInt val : 1;` is one 4-byte word.
        assert_eq!(std::mem::size_of::<CliqueVar>(), 4);
        let v = CliqueVar::new(5, 1);
        assert_eq!((v.col(), v.val(), v.index()), (5, 1, 11));
        assert_eq!(v.word(), 5 | (1 << 31));
        let c = v.complement();
        assert_eq!((c.col(), c.val(), c.word()), (5, 0, 5));
        // bitfield assignment truncation: col keeps 31 bits, val keeps 1 bit
        let t = CliqueVar::new(-1, 3);
        assert_eq!((t.col(), t.val()), (0x7FFF_FFFF, 1));
        let mut s = CliqueVar::new(7, 1);
        s.set_col(0x7FFF_FFFF);
        assert_eq!((s.col(), s.val()), (0x7FFF_FFFF, 1));
        assert_eq!(
            format!("{:?}", CliqueVar::new(3, 0)),
            "HighsCliqueTable_CliqueVar { col: 3, val: 0 }"
        );
    }
}
