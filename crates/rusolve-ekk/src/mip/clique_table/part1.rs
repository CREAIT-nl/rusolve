//! HiGHS source: mip/HighsCliqueTable.cpp (lines 1-1076)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::parallel::highs_combinable::HighsCombinable;
use crate::parallel::highs_parallel::for_each;
use crate::prelude::*;
use crate::util::highs_pdqsort::{pdqsort_branchless_by, pdqsort_by};

/// HiGHS: HighsCliqueTable.cpp::ThreadNeighbourhoodQueryData (lines 475-478)
struct ThreadNeighbourhoodQueryData {
    /// C++: `int64_t numQueries`
    num_queries: i64,
    /// C++: `std::vector<HighsInt> neighbourhoodInds`
    neighbourhood_inds: Vec<HighsInt>,
}

/// HiGHS: HighsCliqueTable.cpp::sortedEdge (lines 27-32)
/// C++: `static std::pair<HighsCliqueTable::CliqueVar, HighsCliqueTable::CliqueVar> sortedEdge(HighsCliqueTable::CliqueVar v1, HighsCliqueTable::CliqueVar v2)`
///
/// PORT NOTE (perf): returns the pair as the two raw 32-bit bitfield words
/// ([`HighsCliqueTable_CliqueVar::word`]), the bytes C++ hashes and compares
/// (`memcmp`) for this trivially copyable key. `sizeTwoCliques` entries are then
/// 12 bytes as in C++ instead of 20, with the same hash and hence the same slot layout.
pub fn sorted_edge(v1: HighsCliqueTable_CliqueVar, v2: HighsCliqueTable_CliqueVar) -> (u32, u32) {
    if v1.col() > v2.col() {
        return (v2.word(), v1.word());
    }

    (v1.word(), v2.word())
}

/// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::CliqueVar::complement applied to every
/// entry of a clique, used by the lambda `complementClique` in addClique (lines 658-662).
fn complement_clique(cliquevars: &mut [HighsCliqueTable_CliqueVar]) {
    for cv in cliquevars.iter_mut() {
        *cv = cv.complement();
    }
}

impl HighsCliqueTable {
    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::unlink (lines 34-43)
    /// C++: `void unlink(HighsInt pos, HighsInt cliqueid)`
    pub fn unlink(&mut self, pos: HighsInt, cliqueid: HighsInt) {
        assert!(pos >= 0);
        let idx = self.cliqueentries[pos as usize].index() as usize;
        self.numcliquesvar[idx] -= 1;

        let cliquelen = self.cliques[cliqueid as usize].end - self.cliques[cliqueid as usize].start;
        if cliquelen == 2 {
            self.inverted_hash_list_size_two[idx].erase(&cliqueid);
        } else {
            self.inverted_hash_list[idx].erase(&cliqueid);
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::link (lines 45-55)
    /// C++: `void link(HighsInt pos, HighsInt cliqueid)`
    pub fn link(&mut self, pos: HighsInt, cliqueid: HighsInt) {
        assert!(pos >= 0);
        debug_assert!(self.col_deleted[self.cliqueentries[pos as usize].col() as usize] == 0);
        let idx = self.cliqueentries[pos as usize].index() as usize;
        self.numcliquesvar[idx] += 1;

        let cliquelen = self.cliques[cliqueid as usize].end - self.cliques[cliqueid as usize].start;
        if cliquelen == 2 {
            self.inverted_hash_list_size_two[idx].insert(HighsHashTableEntry::new(cliqueid, ()));
        } else {
            self.inverted_hash_list[idx].insert(HighsHashTableEntry::new(cliqueid, pos));
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::findCommonCliqueId (lines 57-74)
    /// C++: `HighsInt findCommonCliqueId(int64_t & numQueries, CliqueVar v1, CliqueVar v2) const`
    pub fn find_common_clique_id(
        &self,
        num_queries: &mut i64,
        v1: HighsCliqueTable_CliqueVar,
        v2: HighsCliqueTable_CliqueVar,
    ) -> HighsInt {
        *num_queries += 1;
        if !self.inverted_hash_list_size_two[v1.index() as usize].empty()
            && !self.inverted_hash_list_size_two[v2.index() as usize].empty()
        {
            if let Some(size_two_clique_id) = self.size_two_cliques.find(&sorted_edge(v1, v2)) {
                return *size_two_clique_id;
            }
        }

        let h1 = &self.inverted_hash_list[v1.index() as usize];
        let h2 = &self.inverted_hash_list[v2.index() as usize];

        // HiGHS 1.16 #3327 (h16_hash_common): flag fixed at table construction (see set_h16_flags)
        if let Some(common_clique) = h1.find_common(h2, self.h16_hash_common) {
            return *common_clique.key();
        }
        -1
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::resolveSubstitution (lines 76-81)
    /// C++: `void resolveSubstitution(CliqueVar & v) const`
    pub fn resolve_substitution(&self, v: &mut HighsCliqueTable_CliqueVar) {
        while self.colsubstituted[v.col() as usize] != 0 {
            let subst =
                self.substitutions[(self.colsubstituted[v.col() as usize] - 1) as usize].clone();
            *v = if v.val() == 1 {
                subst.replace
            } else {
                subst.replace.complement()
            };
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::resolveSubstitution (lines 83-93)
    /// C++: `void resolveSubstitution(HighsInt & col, double & val, double & offset) const`
    pub fn resolve_substitution_2(&self, col: &mut HighsInt, val: &mut f64, offset: &mut f64) {
        while self.colsubstituted[*col as usize] != 0 {
            let subst =
                self.substitutions[(self.colsubstituted[*col as usize] - 1) as usize].clone();
            if subst.replace.val() == 0 {
                *offset += *val;
                *val = -*val;
            }
            *col = subst.replace.col() as HighsInt;
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::runCliqueSubsumption (lines 95-129)
    /// C++: `HighsInt runCliqueSubsumption(const HighsDomain & globaldom, std::vector<CliqueVar> & clique)`
    ///
    /// PORT NOTE: C++ `removeCliques` is called from inside `cliqueSubsumption`. The removal is
    /// deferred until that loop has finished and applied in the same order. The loop only reads
    /// the hit counts, origin, numActive and equality of cliques it visits, and each clique is
    /// visited once, so removing a clique earlier or later does not change what is read.
    pub fn run_clique_subsumption(
        &mut self,
        globaldom: &HighsDomain,
        clique: &mut Vec<HighsCliqueTable_CliqueVar>,
    ) -> HighsInt {
        // remove deleted cols
        clique.retain(|clqvar| self.col_deleted[clqvar.col() as usize] == 0);

        // return if clique does not have enough elements
        if clique.len() <= 2 {
            return 0;
        }

        // callback for removing cliques
        let mut nremoved: HighsInt = 0;
        let mut to_remove: Vec<HighsInt> = Vec::new();

        // check clique
        let mut redundant = false;
        let mut dominating_origin: HighsInt = kHighsIInf;
        self.clique_subsumption(
            &clique[..],
            &mut redundant,
            &mut dominating_origin,
            &mut |cliqueid| {
                nremoved += 1;
                to_remove.push(cliqueid);
            },
        );

        for cliqueid in to_remove {
            self.cliques[cliqueid as usize].origin = kHighsIInf;
            self.remove_clique(cliqueid, true);
        }

        if redundant {
            clique.clear();
        }

        if !self.infeasvertexstack.is_empty() {
            clique.retain(|v| !globaldom.is_fixed(v.col() as HighsInt));
        }

        nremoved
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::cliqueSubsumption (lines 131-175)
    /// C++: `void cliqueSubsumption(const std::vector<CliqueVar> & clique, bool & redundant, HighsInt & dominatingOrigin, std::function<void (HighsInt)> removeCliqueCallback)`
    pub fn clique_subsumption(
        &mut self,
        clique: &[HighsCliqueTable_CliqueVar],
        redundant: &mut bool,
        dominating_origin: &mut HighsInt,
        remove_clique_callback: &mut dyn FnMut(HighsInt),
    ) {
        // collect indices of cliques that contain variables from the
        // provided vector
        self.collect_cliques(clique);

        // initialise
        *redundant = false;
        *dominating_origin = kHighsIInf;

        for i in 0..self.cliquehitinds.len() {
            let cliqueid = self.cliquehitinds[i];
            // remember hits and zero out vector
            let hits = self.cliquehits[cliqueid as usize] as HighsInt;
            self.cliquehits[cliqueid as usize] = 0;

            if hits == clique.len() as HighsInt {
                // clique is redundant
                *redundant = true;
                let origin = self.cliques[cliqueid as usize].origin;
                if origin != kHighsIInf && origin != -1 {
                    *dominating_origin = origin;
                }
            } else if self.cliques[cliqueid as usize].num_active() == hits {
                // clique is subset of another clique
                if self.cliques[cliqueid as usize].equality {
                    // subset of clique is an equality clique (stored in the clique table),
                    // thus the remaining variables can be fixed
                    let size_two = self.cliques[cliqueid as usize].end
                        - self.cliques[cliqueid as usize].start
                        == 2;
                    for v in clique {
                        let v_has_clq = if size_two {
                            self.inverted_hash_list_size_two[v.index() as usize].contains(&cliqueid)
                        } else {
                            self.inverted_hash_list[v.index() as usize].contains(&cliqueid)
                        };
                        if !v_has_clq {
                            self.infeasvertexstack.push(*v);
                        }
                    }
                } else {
                    // clique from clique table contains a subset of the variables, thus it
                    // can be removed
                    remove_clique_callback(cliqueid);
                }
            }
        }
        // clear vector of cliques indices
        self.cliquehitinds.clear();
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::collectCliques (lines 177-194)
    /// C++: `void collectCliques(const std::vector<CliqueVar> & clique)`
    ///
    /// PORT NOTE: the C++ `for_each` callbacks only mutate `cliquehits` and `cliquehitinds`, not
    /// the hash tree being iterated. The clique ids are collected first, then applied in the same
    /// order.
    pub fn collect_cliques(&mut self, clique: &[HighsCliqueTable_CliqueVar]) {
        // resize vector
        if self.cliquehits.len() < self.cliques.len() {
            self.cliquehits.resize(self.cliques.len(), 0);
        }

        // collect indices of cliques that contain variables from the
        // provided vector
        for v in clique {
            let mut ids: Vec<HighsInt> = Vec::new();
            self.inverted_hash_list[v.index() as usize].for_each(&mut |cliqueid, _| {
                ids.push(*cliqueid);
            });
            for cliqueid in ids {
                if self.cliquehits[cliqueid as usize] == 0 {
                    self.cliquehitinds.push(cliqueid);
                }
                self.cliquehits[cliqueid as usize] += 1;
            }

            let mut ids_two: Vec<HighsInt> = Vec::new();
            self.inverted_hash_list_size_two[v.index() as usize].for_each(&mut |cliqueid, _| {
                ids_two.push(*cliqueid);
            });
            for cliqueid in ids_two {
                if self.cliquehits[cliqueid as usize] == 0 {
                    self.cliquehitinds.push(cliqueid);
                }
                self.cliquehits[cliqueid as usize] += 1;
            }
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::bronKerboschRecurse (lines 196-304)
    /// C++: `void bronKerboschRecurse(BronKerboschData & data, HighsInt Plen, const CliqueVar * X, HighsInt Xlen) const`
    ///
    /// PORT NOTE: signature changed because the C++ `(const CliqueVar * X, HighsInt Xlen)` pair is
    /// a slice `x`, whose length is `x.len()`.
    pub fn bron_kerbosch_recurse(
        &self,
        data: &mut HighsCliqueTable_BronKerboschData<'_>,
        plen: HighsInt,
        x: &[HighsCliqueTable_CliqueVar],
    ) {
        let mut plen = plen;
        let mut w = data.w_r;
        for i in 0..plen as usize {
            w += data.p[i].weight(data.sol);
        }

        if w < data.min_w - data.feastol {
            return;
        }

        if plen == 0 && x.is_empty() {
            let clique = data.r.clone();

            if data.min_w < w - data.feastol {
                data.maxcliques -= data.cliques.len() as HighsInt;
                data.cliques.clear();
                data.min_w = w;
            }
            data.cliques.push(clique);
            // do not further search for cliques that are violated less than this
            // current clique
            return;
        }

        data.ncalls += 1;

        if data.stop() {
            return;
        }

        let mut pivweight = -1.0;
        let mut pivot = HighsCliqueTable_CliqueVar::from_parts(0, 0);

        for xi in x {
            if xi.weight(data.sol) > pivweight {
                pivweight = xi.weight(data.sol);
                pivot = *xi;
                if pivweight >= 1.0 - data.feastol {
                    break;
                }
            }
        }

        if pivweight < 1.0 - data.feastol {
            for i in 0..plen as usize {
                if data.p[i].weight(data.sol) > pivweight {
                    pivweight = data.p[i].weight(data.sol);
                    pivot = data.p[i];
                    if pivweight >= 1.0 - data.feastol {
                        break;
                    }
                }
            }
        }

        let mut pminus_nu: Vec<HighsCliqueTable_CliqueVar> = Vec::with_capacity(plen as usize);
        self.query_neighbourhood(
            &mut data.neighbourhood_inds,
            &mut data.num_neighbourhood_queries,
            pivot,
            &data.p[..plen as usize],
        );
        data.neighbourhood_inds.push(plen);
        let mut k: HighsInt = 0;
        for i in data.neighbourhood_inds.clone() {
            while k < i {
                pminus_nu.push(data.p[k as usize]);
                k += 1;
            }
            k += 1;
        }

        let sol = data.sol;
        pdqsort_by(&mut pminus_nu, |a, b| {
            (a.weight(sol), a.index()) > (b.weight(sol), b.index())
        });

        let mut local_x: Vec<HighsCliqueTable_CliqueVar> = x.to_vec();

        for v in pminus_nu {
            let new_plen = self.partition_neighbourhood(
                &mut data.neighbourhood_inds,
                &mut data.num_neighbourhood_queries,
                v,
                &mut data.p[..plen as usize],
            );
            let new_xlen = self.partition_neighbourhood(
                &mut data.neighbourhood_inds,
                &mut data.num_neighbourhood_queries,
                v,
                &mut local_x[..],
            );

            // add v to R, update the weight, and do the recursive call
            data.r.push(v);
            let wv = v.weight(data.sol);
            data.w_r += wv;
            self.bron_kerbosch_recurse(data, new_plen, &local_x[..new_xlen as usize]);
            if data.stop() {
                return;
            }

            // remove v from R restore the weight and continue the loop in this call
            data.r.pop();
            data.w_r -= wv;

            w -= wv;
            if w < data.min_w {
                return;
            }
            // find the position of v in the vertices removed from P for the recursive
            // call
            // and also remove it from the set P for this call
            let mut vpos: HighsInt = -1;
            for i in new_plen..plen {
                if data.p[i as usize] == v {
                    vpos = i;
                    break;
                }
            }

            // do the removal by first swapping it to the end of P and reduce the size
            // of P accordingly
            assert!(vpos != -1);

            plen -= 1;
            data.p.swap(vpos as usize, plen as usize);

            local_x.push(v);
        }
    }

    /// HiGHS 1.16 (#3540): HighsCliqueTable.cpp::HighsCliqueTable::addSizeTwoClique (lines 34-51)
    /// C++: `void addSizeTwoClique(HighsInt cliqueid)`
    /// Called only when `h16_clique_links` is set (rule F).
    pub fn add_size_two_clique(&mut self, cliqueid: HighsInt) {
        let start = self.cliques[cliqueid as usize].start;
        if self.cliques[cliqueid as usize].end - start != 2 {
            return;
        }
        let edge = sorted_edge(
            self.cliqueentries[start as usize],
            self.cliqueentries[(start + 1) as usize],
        );
        match self.size_two_cliques.find(&edge).copied() {
            None => {
                // first clique on this edge
                self.size_two_cliques
                    .insert(HighsHashTableEntry::new(edge, cliqueid));
                self.size_two_links[cliqueid as usize] =
                    HighsCliqueTable_SizeTwoLinks { prev: -1, next: -1 };
            }
            Some(prev) => {
                // insert the clique after the stored one
                let next = self.size_two_links[prev as usize].next;
                self.size_two_links[cliqueid as usize] =
                    HighsCliqueTable_SizeTwoLinks { prev, next };
                self.size_two_links[prev as usize].next = cliqueid;
                if next != -1 {
                    self.size_two_links[next as usize].prev = cliqueid;
                }
            }
        }
    }

    /// HiGHS 1.16 (#3540): HighsCliqueTable.cpp::HighsCliqueTable::removeSizeTwoClique (lines 53-72)
    /// C++: `void removeSizeTwoClique(HighsInt cliqueid)`
    /// Called only when `h16_clique_links` is set (rule F).
    pub fn remove_size_two_clique(&mut self, cliqueid: HighsInt) {
        let start = self.cliques[cliqueid as usize].start;
        if self.cliques[cliqueid as usize].end - start != 2 {
            return;
        }
        let prev = self.size_two_links[cliqueid as usize].prev;
        let next = self.size_two_links[cliqueid as usize].next;
        if prev != -1 {
            self.size_two_links[prev as usize].next = next;
        } else {
            // the clique is the stored one, so store the next one if there is one
            let edge = sorted_edge(
                self.cliqueentries[start as usize],
                self.cliqueentries[(start + 1) as usize],
            );
            let size_two_clique_id = self.size_two_cliques.find_mut(&edge);
            // C++ assert (dropped: the oracle is built with NDEBUG)
            debug_assert!(size_two_clique_id.as_deref() == Some(&cliqueid));
            if next != -1 {
                *size_two_clique_id.expect("HiGHS 1.16 #3540: stored size-two clique") = next;
            } else {
                self.size_two_cliques.erase(&edge);
            }
        }
        if next != -1 {
            self.size_two_links[next as usize].prev = prev;
        }
        self.size_two_links[cliqueid as usize] =
            HighsCliqueTable_SizeTwoLinks { prev: -1, next: -1 };
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::doAddClique (lines 334-474)
    /// C++: `void doAddClique(const CliqueVar * cliquevars, HighsInt numcliquevars, bool equality = false, HighsInt origin = kHighsIInf)`
    /// C++ default arguments: `equality = false`, `origin = kHighsIInf`
    ///
    /// PORT NOTE: signature changed because the C++ `(const CliqueVar * cliquevars, HighsInt
    /// numcliquevars)` pair is a slice; its length is `cliquevars.len()`.
    pub fn do_add_clique(
        &mut self,
        cliquevars: &[HighsCliqueTable_CliqueVar],
        equality: bool,
        origin: HighsInt,
    ) {
        // HiGHS 1.16 #3540 (h16_clique_links): read once at the gate
        let h16_clique_links = self.h16_clique_links;
        let numcliquevars = cliquevars.len() as HighsInt;
        let cliqueid: HighsInt;

        if self.freeslots.is_empty() {
            cliqueid = self.cliques.len() as HighsInt;
            self.cliques.push(HighsCliqueTable_Clique::default());
            if h16_clique_links {
                // HiGHS 1.16 #3540 (h16_clique_links): HighsCliqueTable.cpp lines 382 (1.16.0)
                self.size_two_links
                    .push(HighsCliqueTable_SizeTwoLinks { prev: -1, next: -1 });
            }
        } else {
            cliqueid = self.freeslots.pop().unwrap();
        }

        self.cliques[cliqueid as usize].equality = equality;
        self.cliques[cliqueid as usize].origin = origin;

        let max_end: HighsInt;
        let lower = if self.freespaces.is_empty() {
            None
        } else {
            self.freespaces.range((numcliquevars, -1)..).next().copied()
        };
        match lower {
            None => {
                self.cliques[cliqueid as usize].start = self.cliqueentries.len() as HighsInt;
                self.cliques[cliqueid as usize].end =
                    self.cliques[cliqueid as usize].start + numcliquevars;
                max_end = self.cliques[cliqueid as usize].end;
                self.cliqueentries.resize(
                    self.cliques[cliqueid as usize].end as usize,
                    HighsCliqueTable_CliqueVar::default(),
                );
            }
            Some(freespace) => {
                self.freespaces.remove(&freespace);

                self.cliques[cliqueid as usize].start = freespace.1;
                self.cliques[cliqueid as usize].end =
                    self.cliques[cliqueid as usize].start + numcliquevars;
                max_end = self.cliques[cliqueid as usize].start + freespace.0;
            }
        }

        self.cliques[cliqueid as usize].num_zero_fixed = 0;

        let mut fixtozero = false;
        let mut k = self.cliques[cliqueid as usize].start;
        for i in 0..numcliquevars as usize {
            let mut v = cliquevars[i];

            self.resolve_substitution(&mut v);

            if fixtozero {
                self.infeasvertexstack.push(v);
                continue;
            }

            // due to substitutions the variable may occur together with its complement
            // in this clique and we can fix all other variables in the clique to zero:
            //          x + ~x + ... <= 1
            //   <=> x + 1 - x + ... <= 1
            //   <=>             ... <= 0
            let clq_has_vcompl = if numcliquevars == 2 {
                self.inverted_hash_list_size_two[v.complement().index() as usize]
                    .contains(&cliqueid)
            } else {
                self.inverted_hash_list[v.complement().index() as usize].contains(&cliqueid)
            };

            if clq_has_vcompl {
                fixtozero = true;
                for j in self.cliques[cliqueid as usize].start..k {
                    if self.cliqueentries[j as usize].col() != v.col() {
                        self.infeasvertexstack.push(self.cliqueentries[j as usize]);
                    }
                    self.unlink(j, cliqueid);
                }
                k = self.cliques[cliqueid as usize].start;
                continue;
            }

            // due to substitutions the variable may occur twice in this clique and
            // we can fix it to zero:  x + x + ... <= 1  <=>  2x <= 1 <=> x <= 0.5 <=>
            // x = 0
            let inserted = if numcliquevars == 2 {
                self.inverted_hash_list_size_two[v.index() as usize]
                    .insert(HighsHashTableEntry::new(cliqueid, ()))
            } else {
                self.inverted_hash_list[v.index() as usize]
                    .insert(HighsHashTableEntry::new(cliqueid, k))
            };

            if !inserted {
                self.infeasvertexstack.push(v);
                continue;
            }

            self.cliqueentries[k as usize] = v;
            self.numcliquesvar[v.index() as usize] += 1;
            k += 1;
        }

        if max_end > k {
            if self.cliqueentries.len() as HighsInt == max_end {
                self.cliqueentries.truncate(k as usize);
            } else {
                self.freespaces.insert((max_end - k, k));
            }

            if self.cliques[cliqueid as usize].end > k {
                let cstart = self.cliques[cliqueid as usize].start;
                match k - cstart {
                    0 => {
                        // clique empty, so just mark it as deleted
                        self.cliques[cliqueid as usize].start = -1;
                        self.cliques[cliqueid as usize].end = -1;
                        self.freeslots.push(cliqueid);
                        return;
                    }
                    1 => {
                        // size 1 clique is redundant, so unlink the single linked entry
                        // and mark it as deleted
                        self.unlink(cstart, cliqueid);
                        self.cliques[cliqueid as usize].start = -1;
                        self.cliques[cliqueid as usize].end = -1;
                        self.freeslots.push(cliqueid);
                        return;
                    }
                    2 => {
                        // due to substitutions the clique became smaller and is now of size
                        // two as a result we need to link it to the size two cliqueset
                        // instead of the normal cliqueset
                        self.unlink(cstart, cliqueid);
                        self.unlink(cstart + 1, cliqueid);

                        self.cliques[cliqueid as usize].end = k;

                        self.link(cstart, cliqueid);
                        self.link(cstart + 1, cliqueid);
                    }
                    _ => {
                        self.cliques[cliqueid as usize].end = k;
                    }
                }
            }
        }

        let cstart = self.cliques[cliqueid as usize].start;
        let clique_len = self.cliques[cliqueid as usize].end - cstart;
        self.num_entries += clique_len;
        if h16_clique_links {
            // HiGHS 1.16 #3540 (h16_clique_links): HighsCliqueTable.cpp line 510 (1.16.0)
            self.add_size_two_clique(cliqueid);
        } else if clique_len == 2 {
            let key = sorted_edge(
                self.cliqueentries[cstart as usize],
                self.cliqueentries[(cstart + 1) as usize],
            );
            self.size_two_cliques
                .insert(HighsHashTableEntry::new(key, cliqueid));
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::queryNeighbourhood (lines 480-519)
    /// C++: `void queryNeighbourhood(std::vector<HighsInt> & neighbourhoodInds, int64_t & numQueries, CliqueVar v, CliqueVar * q, HighsInt N) const`
    ///
    /// PORT NOTE: signature changed because the C++ `(CliqueVar * q, HighsInt N)` pair is a slice.
    ///
    /// Both C++ branches are ported: the serial loop, and the parallel branch
    /// (`highs::parallel::for_each` grain 10 into a `HighsCombinable`, `combine_each` in
    /// worker-id order, then `pdqsort`).
    pub fn query_neighbourhood(
        &self,
        neighbourhood_inds: &mut Vec<HighsInt>,
        num_queries: &mut i64,
        v: HighsCliqueTable_CliqueVar,
        q: &[HighsCliqueTable_CliqueVar],
    ) {
        neighbourhood_inds.clear();

        if self.num_cliques_2(v) == 0 {
            return;
        }

        // C++ line 482 (`N`): the slice length.
        let n = q.len() as HighsInt;

        // rusolve: same result as the loops below, from the (few) cliques
        // of `v` instead of one hash lookup per candidate.
        if self.query_neighbourhood_by_cliques(neighbourhood_inds, num_queries, v, q) {
            return;
        }

        // C++ lines 487-488: `!allowParallel || numEntries - sizeTwoCliques.size() * 2 <
        // minEntriesForParallelism`; `size()` is `i64` (HighsHash.h 1274), so the arithmetic
        // and the comparison are in i64.
        if !self.allow_parallel
            || (self.num_entries as i64) - self.size_two_cliques.size() * 2
                < self.min_entries_for_parallelism as i64
        {
            // C++ lines 489-491
            for i in 0..n {
                if self.have_common_clique_2(num_queries, v, q[i as usize]) {
                    neighbourhood_inds.push(i);
                }
            }
        } else {
            // C++ lines 493-499: makeHighsCombinable<ThreadNeighbourhoodQueryData>([N]() {...})
            let mut neighbourhood_data = HighsCombinable::with_construct(move || {
                let mut d = ThreadNeighbourhoodQueryData {
                    num_queries: 0,
                    neighbourhood_inds: Vec::new(),
                };
                d.neighbourhood_inds.reserve(n as usize);
                d.num_queries = 0;
                d
            });
            // C++ lines 500-509: highs::parallel::for_each(0, N, [this, &neighbourhoodData, v,
            // q](HighsInt start, HighsInt end) {...}, 10);
            {
                let neighbourhood_data = &neighbourhood_data;
                for_each(
                    0,
                    n,
                    &|start: HighsInt, end: HighsInt| {
                        let d = neighbourhood_data.local();
                        for i in start..end {
                            if self.have_common_clique_2(&mut d.num_queries, v, q[i as usize]) {
                                d.neighbourhood_inds.push(i);
                            }
                        }
                    },
                    10,
                );
            }

            // C++ lines 511-516: combine_each in worker-id order
            neighbourhood_data.combine_each(|d: &mut ThreadNeighbourhoodQueryData| {
                neighbourhood_inds.extend_from_slice(&d.neighbourhood_inds[..]);
                *num_queries += d.num_queries;
            });
            // C++ line 517: pdqsort(begin, end) is pdqsort(begin, end, std::less<HighsInt>()),
            // which takes the branchless loop for an arithmetic type with the default compare
            // (pdqsort.h 495-513).
            pdqsort_branchless_by(&mut neighbourhood_inds[..], |a, b| a < b);
        }
    }

    /// rusolve, not a HiGHS function: `queryNeighbourhood` without one
    /// `findCommonCliqueId` per candidate.
    ///
    /// `haveCommonClique(v, w)` is true when `w` has another column than `v` and some clique holds
    /// both. With `h16_hash_common` (the exact common-key search) that is `w in N(v)`, the union
    /// of the entries of the cliques `v` is linked to, so `N(v)` is collected once from the two
    /// inverted lists and the candidates are filtered against it, in candidate order: same
    /// indices as the scan. `num_queries` grows by one per candidate with another column, as in
    /// the scan. Returns false (nothing written) when the flag is off, or when `v` is in so many
    /// clique entries that the scan is cheaper.
    fn query_neighbourhood_by_cliques(
        &self,
        neighbourhood_inds: &mut Vec<HighsInt>,
        num_queries: &mut i64,
        v: HighsCliqueTable_CliqueVar,
        q: &[HighsCliqueTable_CliqueVar],
    ) -> bool {
        const MAX_NEIGHBOURS: usize = 32;
        if !self.h16_hash_common || q.len() < 64 {
            return false;
        }
        let vcol = v.col();
        let mut nbrs: Vec<u32> = Vec::new();
        let over = std::cell::Cell::new(false);
        let mut visit = |cid: HighsInt| -> bool {
            let cl = &self.cliques[cid as usize];
            for e in &self.cliqueentries[cl.start as usize..cl.end as usize] {
                if e.col() != vcol {
                    if nbrs.len() >= MAX_NEIGHBOURS {
                        over.set(true);
                        return true;
                    }
                    nbrs.push(e.index() as u32);
                }
            }
            false
        };
        let vi = v.index() as usize;
        self.inverted_hash_list_size_two[vi]
            .for_each_until(&mut |cid: &HighsInt, _: &()| visit(*cid));
        if !over.get() {
            self.inverted_hash_list[vi]
                .for_each_until(&mut |cid: &HighsInt, _: &HighsInt| visit(*cid));
        }
        if over.get() {
            return false;
        }
        nbrs.sort_unstable();
        nbrs.dedup();

        let mut queries: i64 = 0;
        for (i, w) in q.iter().enumerate() {
            if w.col() == vcol {
                continue;
            }
            queries += 1;
            if nbrs.binary_search(&(w.index() as u32)).is_ok() {
                neighbourhood_inds.push(i as HighsInt);
            }
        }
        *num_queries += queries;
        true
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::partitionNeighbourhood (lines 521-530)
    /// C++: `HighsInt partitionNeighbourhood(std::vector<HighsInt> & neighbourhoodInds, int64_t & numQueries, CliqueVar v, CliqueVar * q, HighsInt N) const`
    ///
    /// PORT NOTE: signature changed because the C++ `(CliqueVar * q, HighsInt N)` pair is a slice.
    pub fn partition_neighbourhood(
        &self,
        neighbourhood_inds: &mut Vec<HighsInt>,
        num_queries: &mut i64,
        v: HighsCliqueTable_CliqueVar,
        q: &mut [HighsCliqueTable_CliqueVar],
    ) -> HighsInt {
        self.query_neighbourhood(neighbourhood_inds, num_queries, v, q);

        for i in 0..neighbourhood_inds.len() {
            q.swap(i, neighbourhood_inds[i] as usize);
        }

        neighbourhood_inds.len() as HighsInt
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::shrinkToNeighbourhood (lines 532-541)
    /// C++: `HighsInt shrinkToNeighbourhood(std::vector<HighsInt> & neighbourhoodInds, int64_t & numQueries, CliqueVar v, CliqueVar * q, HighsInt N)`
    ///
    /// PORT NOTE: signature changed because the C++ `(CliqueVar * q, HighsInt N)` pair is a slice.
    pub fn shrink_to_neighbourhood(
        &mut self,
        neighbourhood_inds: &mut Vec<HighsInt>,
        num_queries: &mut i64,
        v: HighsCliqueTable_CliqueVar,
        q: &mut [HighsCliqueTable_CliqueVar],
    ) -> HighsInt {
        self.query_neighbourhood(neighbourhood_inds, num_queries, v, q);

        for i in 0..neighbourhood_inds.len() {
            q[i] = q[neighbourhood_inds[i] as usize];
        }

        neighbourhood_inds.len() as HighsInt
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::fixCol (lines 543-554)
    /// C++: `bool fixCol(HighsDomain & globaldom, CliqueVar v, bool doProcessInfeasibleVertices = false)`
    /// C++ default arguments: `do_process_infeasible_vertices = false`
    ///
    /// PORT NOTE: signature changed because C++ `globaldom.fixCol` reaches the solver through the
    /// domain, so the solver is passed alongside the domain (M1 context).
    pub fn fix_col(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldom: &mut HighsDomain,
        v: HighsCliqueTable_CliqueVar,
        do_process_infeasible_vertices: bool,
    ) -> bool {
        let wasfixed = globaldom.is_fixed(v.col() as HighsInt);
        lend!(
            ms.mipdata_.cliquetable,
            self,
            globaldom.fix_col(
                ms,
                v.col() as HighsInt,
                (1 - v.val() as HighsInt) as f64,
                HighsDomain_Reason::unspecified(),
            )
        );
        if globaldom.infeasible() {
            return false;
        }
        if !wasfixed {
            self.nfixings += 1;
            self.infeasvertexstack.push(v);
            if do_process_infeasible_vertices {
                self.process_infeasible_vertices(ms, globaldom);
            }
        }
        true
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::processNewEdge (lines 556-646)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::processNewEdge (lines 556-602)
    /// C++: `bool processNewEdge(HighsDomain & globaldom, CliqueVar v1, CliqueVar v2)`
    ///
    /// PORT NOTE: signature changed because `fixCol` and `foundCover` reach the solver through the
    /// domain (M1 context).
    pub fn process_new_edge(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldom: &mut HighsDomain,
        v1: HighsCliqueTable_CliqueVar,
        v2: HighsCliqueTable_CliqueVar,
    ) -> bool {
        let mut v1 = v1;
        let mut v2 = v2;

        if v1.col() == v2.col() {
            if v1.val() == v2.val() {
                self.fix_col(ms, globaldom, v1, true);
                return false;
            }

            return true;
        }

        if self.have_common_clique(v1.complement(), v2) {
            self.fix_col(ms, globaldom, v2, true);
            return false;
        } else if self.have_common_clique(v2.complement(), v1) {
            self.fix_col(ms, globaldom, v1, true);
            return false;
        }

        // return if complements are not in a clique
        if !self.found_cover(ms, globaldom, v1.complement(), v2.complement()) {
            return false;
        }
        if globaldom.infeasible() {
            return true;
        }

        self.found_cover(ms, globaldom, v1, v2);
        if globaldom.is_fixed(v1.col() as HighsInt)
            || globaldom.is_fixed(v2.col() as HighsInt)
            || globaldom.infeasible()
        {
            return true;
        }

        let substitution = if v2.col() < v1.col() {
            if v1.val() == 1 {
                v2 = v2.complement();
            }

            HighsCliqueTable_Substitution {
                substcol: v1.col() as HighsInt,
                replace: v2,
            }
        } else {
            if v2.val() == 1 {
                v1 = v1.complement();
            }

            HighsCliqueTable_Substitution {
                substcol: v2.col() as HighsInt,
                replace: v1,
            }
        };

        // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp line 598 at 7b890c255
        if self.x1b_presolve_clique_table {
            self.record_substitution(substitution);
        } else {
            self.substitutions.push(substitution.clone());
            self.colsubstituted[substitution.substcol as usize] =
                self.substitutions.len() as HighsInt;

            self.replace_substituted_lists(
                HighsCliqueTable_CliqueVar::new(substitution.substcol, 1),
                substitution.replace,
            );
            self.replace_substituted_lists(
                HighsCliqueTable_CliqueVar::new(substitution.substcol, 0),
                substitution.replace.complement(),
            );
        }

        true
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::processNewEdge, lambda `replace` (lines 598-636)
    ///
    /// PORT NOTE: the C++ lambda is a method. Each `for_each` over the substituted variable's list
    /// is collected first, then applied. The C++ callbacks write only to the replacement list,
    /// `cliqueentries` and the size-two hash table, and never to the list being iterated, so the
    /// order of effects is the same.
    fn replace_substituted_lists(
        &mut self,
        substituted_var: HighsCliqueTable_CliqueVar,
        replacement_var: HighsCliqueTable_CliqueVar,
    ) {
        let subst_idx = substituted_var.index() as usize;
        let repl_idx = replacement_var.index() as usize;

        self.numcliquesvar[repl_idx] += self.numcliquesvar[subst_idx];
        self.numcliquesvar[subst_idx] = 0;

        let mut moved: Vec<(HighsInt, HighsInt)> = Vec::new();
        self.inverted_hash_list[subst_idx].for_each(&mut |cliqueid, location| {
            moved.push((*cliqueid, *location));
        });
        for (cliqueid, location) in moved {
            self.inverted_hash_list[repl_idx].insert(HighsHashTableEntry::new(cliqueid, location));
            self.cliqueentries[location as usize] = replacement_var;
        }

        self.inverted_hash_list[subst_idx].clear();

        let mut moved_two: Vec<HighsInt> = Vec::new();
        self.inverted_hash_list_size_two[subst_idx].for_each(&mut |cliqueid, _| {
            moved_two.push(*cliqueid);
        });
        for cliqueid in moved_two {
            let mut pos = self.cliques[cliqueid as usize].start;
            let mut other_pos = pos + 1;

            if self.cliqueentries[other_pos as usize] == substituted_var {
                std::mem::swap(&mut pos, &mut other_pos);
            }

            self.inverted_hash_list_size_two[repl_idx]
                .insert(HighsHashTableEntry::new(cliqueid, ()));
            // HiGHS 1.16 #3540 (h16_clique_links): read once at the gate
            let h16_clique_links = self.h16_clique_links;
            if h16_clique_links {
                // HiGHS 1.16 #3540 (h16_clique_links): HighsCliqueTable.cpp lines 670-672 (1.16.0)
                self.remove_size_two_clique(cliqueid);
                self.cliqueentries[pos as usize] = replacement_var;
                self.add_size_two_clique(cliqueid);
            } else {
                self.cliqueentries[pos as usize] = replacement_var;

                let other = self.cliqueentries[other_pos as usize];
                self.size_two_cliques
                    .erase(&sorted_edge(substituted_var, other));
                self.size_two_cliques.insert(HighsHashTableEntry::new(
                    sorted_edge(replacement_var, other),
                    cliqueid,
                ));
            }
        }

        self.inverted_hash_list_size_two[subst_idx].clear();
    }

    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::replaceLiteral (lines 604-641)
    /// C++: `void replaceLiteral(CliqueVar substitutedVar, CliqueVar replacementVar)`
    ///
    /// #3312 extracts the 1.15.1 `processNewEdge` lambda `replace` unchanged, so this is the
    /// same body ([`HighsCliqueTable::replace_substituted_lists`], which also carries the
    /// 1.16.0 #3540 size-two-link form). Called by `HPresolveCliqueTable::substituteCol`.
    pub fn replace_literal(
        &mut self,
        substituted_var: HighsCliqueTable_CliqueVar,
        replacement_var: HighsCliqueTable_CliqueVar,
    ) {
        self.replace_substituted_lists(substituted_var, replacement_var);
    }

    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::recordSubstitution (lines 643-651)
    /// C++: `void recordSubstitution(Substitution substitution)`
    pub fn record_substitution(&mut self, substitution: HighsCliqueTable_Substitution) {
        self.substitutions.push(substitution.clone());
        self.colsubstituted[substitution.substcol as usize] = self.substitutions.len() as HighsInt;

        self.replace_literal(
            HighsCliqueTable_CliqueVar::new(substitution.substcol, 1),
            substitution.replace,
        );
        self.replace_literal(
            HighsCliqueTable_CliqueVar::new(substitution.substcol, 0),
            substitution.replace.complement(),
        );
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::addClique (lines 648-812)
    /// C++: `void addClique(const HighsMipSolver & mipsolver, CliqueVar * cliquevars, HighsInt numcliquevars, bool equality = false, HighsInt origin = kHighsIInf)`
    /// C++ default arguments: `equality = false`, `origin = kHighsIInf`
    ///
    /// PORT NOTE: signature changed because C++ addClique reads and writes the global domain,
    /// the node queue and the pruned tree weight, so the solver is taken mutably (M1). The domain
    /// is detached for the duration of the body. The `(CliqueVar * cliquevars, numcliquevars)`
    /// pair is a slice.
    pub fn add_clique(
        &mut self,
        ms: &mut HighsMipSolver,
        cliquevars: &mut [HighsCliqueTable_CliqueVar],
        equality: bool,
        origin: HighsInt,
    ) {
        ms.mipdata_.debug_solution.check_clique(&cliquevars[..]);
        detach!(ms.mipdata_.domains[0], |globaldom| self
            .add_clique_inner(ms, globaldom, cliquevars, equality, origin));
    }

    /// PORT NOTE: body of [`HighsCliqueTable::add_clique`] with the global domain passed in
    /// detached. The C++ early `return`s are kept here so that the domain is always re-attached
    /// by the caller.
    fn add_clique_inner(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldom: &mut HighsDomain,
        cliquevars: &mut [HighsCliqueTable_CliqueVar],
        equality: bool,
        origin: HighsInt,
    ) {
        // resolve substitutions
        for i in 0..cliquevars.len() {
            self.resolve_substitution(&mut cliquevars[i]);
        }

        // initialise flag
        let mut has_new_edge = false;

        // check if all variables can be fixed or infeasibility is detected
        if self.check_clique_complemented(ms, globaldom, cliquevars, &mut has_new_edge, false) {
            return;
        }
        if globaldom.infeasible() {
            return;
        }

        if cliquevars.len() == 2 && equality {
            // try complemented clique
            if self.check_clique_complemented(ms, globaldom, cliquevars, &mut has_new_edge, true) {
                return;
            }
            if globaldom.infeasible() {
                return;
            }
        }

        if !has_new_edge && origin == kHighsIInf {
            return;
        }

        // PORT NOTE: std::remove_if over the fixed columns, written as an in-place compaction.
        // Only the first `numunfixed` entries are read afterwards, as in C++.
        let mut numunfixed = 0usize;
        for i in 0..cliquevars.len() {
            if !globaldom.is_fixed(cliquevars[i].col() as HighsInt) {
                cliquevars[numunfixed] = cliquevars[i];
                numunfixed += 1;
            }
        }
        if numunfixed < 2 {
            return;
        }

        self.do_add_clique(&cliquevars[..numunfixed], equality, origin);
        self.process_infeasible_vertices(ms, globaldom);
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::addClique, lambda `checkClique` (lines 743-751)
    fn check_clique_complemented(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldom: &mut HighsDomain,
        cliquevars: &mut [HighsCliqueTable_CliqueVar],
        has_new_edge: &mut bool,
        complement: bool,
    ) -> bool {
        if complement {
            complement_clique(cliquevars);
        }
        let done = self.fix_all_vars_in_clique(ms, globaldom, cliquevars, has_new_edge);
        if complement {
            complement_clique(cliquevars);
        }
        done
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::addClique, lambda `fixAllVarsInClique` (lines 667-734)
    fn fix_all_vars_in_clique(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldom: &mut HighsDomain,
        cliquevars: &mut [HighsCliqueTable_CliqueVar],
        has_new_edge: &mut bool,
    ) -> bool {
        let numcliquevars = cliquevars.len();

        for i in 0..numcliquevars {
            let cvi = cliquevars[i];
            if !globaldom.is_fixed(cvi.col() as HighsInt)
                || (cvi.val() as f64) != globaldom.col_lower_[cvi.col() as usize]
            {
                continue;
            }
            // column is fixed to 1, every other entry can be fixed to zero
            for k in 0..numcliquevars {
                if k == i {
                    continue;
                }
                if !self.fix_col(ms, globaldom, cliquevars[k], false) {
                    return false;
                }
            }
            self.process_infeasible_vertices(ms, globaldom);
            return true;
        }

        if numcliquevars > 100 {
            return false;
        }

        // todo, sort new clique to allow log n lookup of membership in size by
        // binary search

        for i in 0..numcliquevars.saturating_sub(1) {
            if globaldom.is_fixed(cliquevars[i].col() as HighsInt) {
                continue;
            }
            if self.num_cliques_2(cliquevars[i]) == 0
                && self.num_cliques_2(cliquevars[i].complement()) == 0
            {
                *has_new_edge = true;
                continue;
            }

            for j in i + 1..numcliquevars {
                if globaldom.is_fixed(cliquevars[j].col() as HighsInt) {
                    continue;
                }

                if self.have_common_clique(cliquevars[i], cliquevars[j]) {
                    continue;
                }
                // todo: Instead of haveCommonClique use findCommonClique. If the
                // common clique is smaller than this clique check if it is a subset
                // of this clique. If it is a subset remove the clique and iterate the
                // process until either a common clique that is not a subset of this
                // one is found, or no common clique exists anymore in which case we
                // proceed with the code below and set hasNewEdge to true

                *has_new_edge = true;

                let iscover = self.process_new_edge(ms, globaldom, cliquevars[i], cliquevars[j]);
                if globaldom.infeasible() {
                    return false;
                }

                if !ms.mipdata_.nodequeue.empty() {
                    prune_overlapping_nodes(ms, cliquevars[i], cliquevars[j]);
                }

                if iscover {
                    for k in 0..numcliquevars {
                        if k == i || k == j {
                            continue;
                        }
                        if !self.fix_col(ms, globaldom, cliquevars[k], false) {
                            return false;
                        }
                    }
                    self.process_infeasible_vertices(ms, globaldom);
                    return true;
                }
            }
        }
        false
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::removeClique (lines 814-837)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::removeClique (lines 818-842)
    /// C++: `void removeClique(HighsInt cliqueid, bool recordDeletedRow = true)`
    /// C++ default arguments: `recordDeletedRow = true` (the 1.15.1 behaviour; only
    /// `mip_presolve_clique_table` code passes `false`).
    pub fn remove_clique(&mut self, cliqueid: HighsInt, record_deleted_row: bool) {
        if record_deleted_row {
            if self.cliques[cliqueid as usize].origin != kHighsIInf
                && self.cliques[cliqueid as usize].origin != -1
            {
                self.deletedrows
                    .push(self.cliques[cliqueid as usize].origin);
            }
        }

        let start = self.cliques[cliqueid as usize].start;
        assert!(start != -1);
        let end = self.cliques[cliqueid as usize].end;
        let len = end - start;
        // HiGHS 1.16 #3540 (h16_clique_links): read once at the gate
        let h16_clique_links = self.h16_clique_links;
        if h16_clique_links {
            // HiGHS 1.16 #3540 (h16_clique_links): HighsCliqueTable.cpp line 862 (1.16.0)
            self.remove_size_two_clique(cliqueid);
        } else if len == 2 {
            let key = sorted_edge(
                self.cliqueentries[start as usize],
                self.cliqueentries[(start + 1) as usize],
            );
            self.size_two_cliques.erase(&key);
        }

        for i in start..end {
            self.unlink(i, cliqueid);
        }

        self.freeslots.push(cliqueid);
        self.freespaces.insert((len, start));

        self.cliques[cliqueid as usize].start = -1;
        self.cliques[cliqueid as usize].end = -1;
        self.num_entries -= len;
    }

    /// HiGHS latest (#3220, b7d4e98c5): HighsCliqueTable.cpp::HighsCliqueTable::fixLastActiveAndRemove (lines 839-849)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::fixLastActiveAndRemove (lines 844-859)
    pub fn fix_last_active_and_remove(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldom: &mut HighsDomain,
        cliqueid: HighsInt,
    ) {
        // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 844-859 at 7b890c255
        let clqp = self.x1b_presolve_clique_table;
        if self.cliques[cliqueid as usize].equality
            && self.cliques[cliqueid as usize].num_active() == 1
        {
            let start = self.cliques[cliqueid as usize].start;
            let end = self.cliques[cliqueid as usize].end;
            if clqp {
                for i in start..end {
                    let ce = self.cliqueentries[i as usize];
                    if self.col_deleted[ce.col() as usize] != 0 {
                        continue;
                    }
                    if !globaldom.is_fixed(ce.col() as HighsInt) {
                        self.fix_col(ms, globaldom, ce.complement(), false);
                    } else if globaldom.col_lower_[ce.col() as usize] != ce.val() as f64 {
                        globaldom.mark_infeasible(HighsDomain_Reason::unspecified());
                    }
                    break;
                }
            } else {
                for i in start..end {
                    if !globaldom.is_fixed(self.cliqueentries[i as usize].col() as HighsInt) {
                        let v = self.cliqueentries[i as usize].complement();
                        self.fix_col(ms, globaldom, v, false);
                        break;
                    }
                }
            }
        }
        self.remove_clique(cliqueid, true);
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::extractCliques (lines 839-977)
    /// C++: `void extractCliques(const HighsMipSolver & mipsolver, std::vector<HighsInt> & inds, std::vector<double> & vals, std::vector<int8_t> & complementation, double rhs, HighsInt nbin, std::vector<HighsInt> & perm, std::vector<CliqueVar> & clique, double feastol)`
    ///
    /// PORT NOTE: signature changed because C++ extractCliques adds variable bounds to the
    /// implications and adds cliques through the solver, so the solver is taken mutably (M1).
    /// The C++ `nbin` argument is overwritten before it is read, so its value is ignored here
    /// (`_nbin`).
    pub fn extract_cliques(
        &mut self,
        ms: &mut HighsMipSolver,
        inds: &mut Vec<HighsInt>,
        vals: &mut Vec<f64>,
        complementation: &mut Vec<i8>,
        rhs: f64,
        _nbin: HighsInt,
        perm: &mut Vec<HighsInt>,
        clique: &mut Vec<HighsCliqueTable_CliqueVar>,
        feastol: f64,
    ) {
        perm.clear();
        perm.resize(inds.len(), 0);
        for (i, p) in perm.iter_mut().enumerate() {
            *p = i as HighsInt;
        }

        let binaryend = cpp_partition(&mut perm[..], |&pos| {
            ms.mipdata_.domains[0].is_binary(&*ms, inds[pos as usize])
        });
        let mut nbin = binaryend as HighsInt;
        let ntotal = perm.len() as HighsInt;

        // if not all variables are binary, we extract variable upper and lower bounds
        // constraints on the non-binary variable for each binary variable in the
        // constraint
        if nbin < ntotal {
            for i in 0..nbin as usize {
                let bincol = inds[perm[i] as usize];
                let impliedub = HighsCDouble::from(rhs) - vals[perm[i] as usize];
                if ms.mipdata_.implications.too_many_var_bounds() {
                    break;
                }
                for j in nbin as usize..ntotal as usize {
                    let col = inds[perm[j] as usize];
                    if ms.mipdata_.domains[0].is_fixed(col) {
                        continue;
                    }

                    let colub = HighsCDouble::from(ms.mipdata_.domains[0].col_upper_[col as usize])
                        - ms.mipdata_.domains[0].col_lower_[col as usize];
                    let mut implcolub = impliedub / vals[perm[j] as usize];
                    if ms.is_col_integral(col) {
                        let t = f64::from(implcolub) + ms.mipdata_.feastol;
                        implcolub = HighsCDouble::from(t.floor());
                    }

                    if implcolub < colub - feastol {
                        let mut coef: HighsCDouble;
                        let mut constant: HighsCDouble;

                        if complementation[perm[i] as usize] == -1 {
                            coef = colub - implcolub;
                            constant = implcolub;
                        } else {
                            coef = implcolub - colub;
                            constant = colub;
                        }

                        if complementation[perm[j] as usize] == -1 {
                            constant -= ms.mipdata_.domains[0].col_upper_[col as usize];
                            let neg_coef = -f64::from(coef);
                            let neg_constant = -f64::from(constant);
                            detach!(ms.mipdata_.implications, |implics| implics.add_vlb(
                                ms,
                                col,
                                bincol,
                                neg_coef,
                                neg_constant
                            ));
                        } else {
                            constant += ms.mipdata_.domains[0].col_lower_[col as usize];
                            let pos_coef = f64::from(coef);
                            let pos_constant = f64::from(constant);
                            detach!(ms.mipdata_.implications, |implics| implics.add_vub(
                                ms,
                                col,
                                bincol,
                                pos_coef,
                                pos_constant
                            ));
                        }
                    }
                }
            }
        }

        // only one binary means we do have no cliques
        if nbin <= 1 {
            return;
        }

        pdqsort_by(&mut perm[..nbin as usize], |&p1, &p2| {
            (vals[p1 as usize], p1) > (vals[p2 as usize], p2)
        });
        // check if any cliques exists
        if vals[perm[0] as usize] + vals[perm[1] as usize] <= rhs + feastol {
            return;
        }

        // check if this is a set packing constraint (or easily transformable
        // into one)
        let last_bin = vals[perm[(nbin - 1) as usize] as usize];
        if (vals[perm[0] as usize] - last_bin).abs() <= feastol && rhs < 2.0 * last_bin - feastol {
            // the coefficients on the binary variables are all equal and the
            // right hand side is strictly below two times the coefficient value.
            // Therefore the constraint can be transformed into a set packing
            // constraint by relaxing out all non-binary variables (if any),
            // dividing by the coefficient value of the binary variables, and then
            // possibly rounding down the right hand side
            clique.clear();

            for j in 0..nbin as usize {
                let pos = perm[j] as usize;
                if complementation[pos] == -1 {
                    clique.push(HighsCliqueTable_CliqueVar::new(inds[pos], 0));
                } else {
                    clique.push(HighsCliqueTable_CliqueVar::new(inds[pos], 1));
                }
            }

            let clique_len = nbin as usize;
            self.add_clique(ms, &mut clique[..clique_len], false, kHighsIInf);
            if ms.mipdata_.domains[0].infeasible() {
                return;
            }
            return;
        }

        for k in (1..nbin as usize).rev() {
            let mincliqueval = rhs - vals[perm[k] as usize] + feastol;
            let cliqueend = perm[..k].partition_point(|&p| vals[p as usize] > mincliqueval);

            // no clique for this variable
            if cliqueend == 0 {
                continue;
            }

            clique.clear();

            for j in 0..cliqueend {
                let pos = perm[j] as usize;
                if complementation[pos] == -1 {
                    clique.push(HighsCliqueTable_CliqueVar::new(inds[pos], 0));
                } else {
                    clique.push(HighsCliqueTable_CliqueVar::new(inds[pos], 1));
                }
            }

            if complementation[perm[k] as usize] == -1 {
                clique.push(HighsCliqueTable_CliqueVar::new(inds[perm[k] as usize], 0));
            } else {
                clique.push(HighsCliqueTable_CliqueVar::new(inds[perm[k] as usize], 1));
            }

            // further cliques are just subsets of this clique

            if clique.len() >= 2 {
                let clique_len = clique.len();
                self.add_clique(ms, &mut clique[..clique_len], false, kHighsIInf);
                if ms.mipdata_.domains[0].infeasible() {
                    return;
                }
            }

            // further cliques are just subsets of this clique
            if cliqueend == k {
                return;
            }
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::cliquePartition (lines 979-1006)
    /// C++: `void cliquePartition(std::vector<CliqueVar> & clqVars, std::vector<HighsInt> & partitionStart)`
    ///
    /// PORT NOTE: `numNeighbourhoodQueries` is a member passed by reference into a const method, so
    /// its value is copied in before each query and copied back afterwards.
    pub fn clique_partition(
        &mut self,
        clq_vars: &mut Vec<HighsCliqueTable_CliqueVar>,
        partition_start: &mut Vec<HighsInt>,
    ) {
        let n = clq_vars.len() as HighsInt;
        self.randgen.shuffle(&mut clq_vars[..], n);

        let mut neighbourhood_inds: Vec<HighsInt> = Vec::with_capacity(clq_vars.len());

        let num_clq_vars = clq_vars.len() as HighsInt;
        partition_start.clear();
        partition_start.reserve(clq_vars.len());
        let mut extension_end = num_clq_vars;
        partition_start.push(0);
        for i in 0..num_clq_vars {
            if i == extension_end {
                partition_start.push(i);
                extension_end = num_clq_vars;
            }
            let v = clq_vars[i as usize];
            let extension_start = i + 1;
            let mut num_queries = self.num_neighbourhood_queries;
            let count = self.partition_neighbourhood(
                &mut neighbourhood_inds,
                &mut num_queries,
                v,
                &mut clq_vars[extension_start as usize..extension_end as usize],
            );
            self.num_neighbourhood_queries = num_queries;
            extension_end = count + extension_start;
        }

        partition_start.push(num_clq_vars);
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::cliquePartition (lines 1008-1054)
    /// C++: `void cliquePartition(const std::vector<double> & objective, std::vector<CliqueVar> & clqVars, std::vector<HighsInt> & partitionStart)`
    ///
    /// PORT NOTE: `numNeighbourhoodQueries` is copied in and out around each const query, as in
    /// [`HighsCliqueTable::clique_partition`]. The objective factor `(2 * val - 1)` is computed as
    /// a signed integer, since C++ promotes the 1-bit `val` bitfield to `int`.
    pub fn clique_partition_2(
        &mut self,
        objective: &[f64],
        clq_vars: &mut Vec<HighsCliqueTable_CliqueVar>,
        partition_start: &mut Vec<HighsInt>,
    ) {
        let n = clq_vars.len() as HighsInt;
        self.randgen.shuffle(&mut clq_vars[..], n);

        let cmp = |v1: &HighsCliqueTable_CliqueVar, v2: &HighsCliqueTable_CliqueVar| {
            (2 * v1.val() as i64 - 1) as f64 * objective[v1.col() as usize]
                > (2 * v2.val() as i64 - 1) as f64 * objective[v2.col() as usize]
        };
        pdqsort_branchless_by(&mut clq_vars[..], cmp);

        let mut neighbourhood_inds: Vec<HighsInt> = Vec::with_capacity(clq_vars.len());

        let num_clq_vars = clq_vars.len() as HighsInt;
        partition_start.clear();
        partition_start.reserve(clq_vars.len());
        let mut extension_end = num_clq_vars;
        partition_start.push(0);
        let mut last_swapped_index: HighsInt = 0;
        for i in 0..num_clq_vars {
            if i == extension_end {
                partition_start.push(i);
                extension_end = num_clq_vars;
                if last_swapped_index >= i {
                    pdqsort_branchless_by(
                        &mut clq_vars[i as usize..(last_swapped_index + 1) as usize],
                        cmp,
                    );
                }
                last_swapped_index = 0;
            }
            let v = clq_vars[i as usize];
            let extension_start = i + 1;
            let mut num_queries = self.num_neighbourhood_queries;
            let count = self.partition_neighbourhood(
                &mut neighbourhood_inds,
                &mut num_queries,
                v,
                &mut clq_vars[extension_start as usize..extension_end as usize],
            );
            self.num_neighbourhood_queries = num_queries;
            extension_end = count + extension_start;
            if let Some(&back) = neighbourhood_inds.last() {
                last_swapped_index = cpp_max(back + extension_start, last_swapped_index);
            }
        }

        partition_start.push(num_clq_vars);
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsCliqueTable.cpp::HighsCliqueTable::extendClique (lines 988-1005)
    /// C++: `template <bool Sort, typename Comparator> HighsInt extendClique(std::vector<HighsInt>& neighbourhoodInds, std::vector<CliqueVar>& clqVars, HighsInt seedPos, HighsInt extensionEnd, const Comparator& candidateOrder)`
    ///
    /// PORT NOTE: only the sorting overload is ported (it is the one `cliqueCover` uses). The C++
    /// `numNeighbourhoodQueries` member is passed as `num_queries`, a copy the caller writes back.
    pub fn extend_clique<
        F: FnMut(&HighsCliqueTable_CliqueVar, &HighsCliqueTable_CliqueVar) -> bool,
    >(
        &self,
        neighbourhood_inds: &mut Vec<HighsInt>,
        num_queries: &mut i64,
        clq_vars: &mut [HighsCliqueTable_CliqueVar],
        seed_pos: HighsInt,
        mut extension_end: HighsInt,
        candidate_order: &mut F,
    ) -> HighsInt {
        let mut i = seed_pos;
        while i < extension_end {
            pdqsort_branchless_by(
                &mut clq_vars[i as usize..extension_end as usize],
                &mut *candidate_order,
            );
            let v = clq_vars[i as usize];
            let extension_start = i + 1;
            extension_end = self.partition_neighbourhood(
                neighbourhood_inds,
                num_queries,
                v,
                &mut clq_vars[extension_start as usize..extension_end as usize],
            ) + extension_start;
            i += 1;
        }
        extension_end
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsCliqueTable.cpp::HighsCliqueTable::cliqueCover (lines 1082-1134)
    /// C++: `void cliqueCover(std::vector<CliqueVar>& clqVars, std::vector<std::vector<CliqueVar>>& cover)`
    pub fn clique_cover(
        &mut self,
        clq_vars: &mut Vec<HighsCliqueTable_CliqueVar>,
        cover: &mut Vec<Vec<HighsCliqueTable_CliqueVar>>,
    ) {
        cover.clear();
        if clq_vars.is_empty() {
            return;
        }

        let mut neighbourhood_inds: Vec<HighsInt> = Vec::with_capacity(clq_vars.len());

        let num_clq_vars = clq_vars.len() as HighsInt;
        let mut pos: HighsInt = 0;

        while pos < num_clq_vars {
            let mut candidate_order =
                |a: &HighsCliqueTable_CliqueVar, b: &HighsCliqueTable_CliqueVar| {
                    self.numcliquesvar[a.index() as usize] > self.numcliquesvar[b.index() as usize]
                };

            let mut num_queries = self.num_neighbourhood_queries;
            let clique_end = self.extend_clique(
                &mut neighbourhood_inds,
                &mut num_queries,
                &mut clq_vars[..],
                pos,
                num_clq_vars,
                &mut candidate_order,
            );
            self.num_neighbourhood_queries = num_queries;

            let mut clique: Vec<HighsCliqueTable_CliqueVar> =
                clq_vars[pos as usize..clique_end as usize].to_vec();

            for i in 0..pos as usize {
                let v = clq_vars[i];
                let all_common = clique
                    .iter()
                    .all(|clqvar| self.have_common_clique(v, *clqvar));
                if all_common {
                    clique.push(v);
                }
            }

            cover.push(clique);
            pos = clique_end;
        }

        // count variable occurrences
        let mut count: HighsHashTable<HighsInt, HighsInt> = HighsHashTable::new();
        for clique in cover.iter() {
            for var in clique.iter() {
                *count.index_or_insert(&var.index()) += 1;
            }
        }

        // remove cliques where every member appears in at least one other clique
        let mut kept: usize = 0;
        for idx in 0..cover.len() {
            let all_redundant = cover[idx]
                .iter()
                .all(|var| *count.index_or_insert(&var.index()) >= 2);
            if all_redundant {
                for var in cover[idx].iter() {
                    *count.index_or_insert(&var.index()) -= 1;
                }
            } else {
                cover.swap(kept, idx);
                kept += 1;
            }
        }
        cover.truncate(kept);
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::foundCover (lines 1056-1076)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::foundCover (lines 1148-1170)
    /// C++: `bool foundCover(HighsDomain & globaldom, CliqueVar v1, CliqueVar v2)`
    ///
    /// PORT NOTE: signature changed because `fixCol` reaches the solver through the domain (M1).
    pub fn found_cover(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldom: &mut HighsDomain,
        v1: HighsCliqueTable_CliqueVar,
        v2: HighsCliqueTable_CliqueVar,
    ) -> bool {
        // HiGHS #3312 (mip_presolve_clique_table): read once at the gate
        let clqp = self.x1b_presolve_clique_table;
        let mut commonclique = self.find_common_clique_id_2(v1, v2);
        if commonclique == -1 {
            return false;
        }

        while commonclique != -1 {
            let start = self.cliques[commonclique as usize].start;
            let end = self.cliques[commonclique as usize].end;

            for i in start..end {
                let entry = self.cliqueentries[i as usize];
                if clqp {
                    // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 1158-1160 at 7b890c255
                    if entry == v1 || entry == v2 || self.col_deleted[entry.col() as usize] != 0 {
                        continue;
                    }
                } else if entry == v1 || entry == v2 {
                    continue;
                }
                if !self.fix_col(ms, globaldom, entry, false) {
                    return true;
                }
            }

            self.remove_clique(commonclique, true);
            commonclique = self.find_common_clique_id_2(v1, v2);
        }

        self.process_infeasible_vertices(ms, globaldom);
        true
    }
}

/// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::addClique, the node-pruning block (lines 700-727)
///
/// PORT NOTE: the C++ takes `lower_bound` and `upper_bound` of `(val, kHighsIInf)` in each node set,
/// which is the range of entries equal to that key. The entries are collected as a snapshot before
/// any pruning. `pruneNode` unlinks only the node it is given, and the walk has already advanced past
/// it, so the snapshot gives the same result as the C++ iterators.
fn prune_overlapping_nodes(
    ms: &mut HighsMipSolver,
    cv_i: HighsCliqueTable_CliqueVar,
    cv_j: HighsCliqueTable_CliqueVar,
) {
    let (a, b) = {
        let nq = &ms.mipdata_.nodequeue;
        let v1_set = if cv_i.val() == 1 {
            nq.get_up_nodes(cv_i.col() as HighsInt)
        } else {
            nq.get_down_nodes(cv_i.col() as HighsInt)
        };
        let v2_set = if cv_j.val() == 1 {
            nq.get_up_nodes(cv_j.col() as HighsInt)
        } else {
            nq.get_down_nodes(cv_j.col() as HighsInt)
        };
        if v1_set.is_empty() || v2_set.is_empty() {
            return;
        }
        // HiGHS quirk: C++ (HighsCliqueTable.cpp:730-737) takes both `lower_bound` and
        // `upper_bound` at the same key `(val, kHighsIInf)`. The range [lower, upper) therefore
        // holds only an element equal to `(val, kHighsIInf)`. Node ids stay below kHighsIInf, so
        // the range is in practice always empty and the pruning below never runs. The comment
        // above the C++ code suggests that `lower_bound(val, -inf)` was intended. We reproduce
        // the code as written: `range(key..=key)` is exactly [lower_bound, upper_bound).
        let key1 = (OrderedF64(cv_i.val() as f64), kHighsIInf as i64);
        let key2 = (OrderedF64(cv_j.val() as f64), kHighsIInf as i64);
        let a: Vec<i64> = v1_set.range(key1..=key1).map(|e| e.1).collect();
        let b: Vec<i64> = v2_set.range(key2..=key2).map(|e| e.1).collect();
        (a, b)
    };

    if a.is_empty() || b.is_empty() {
        return;
    }
    // itV1 != endV1 && itV2 != endV2 && (itV1->second <= prev(endV2)->second || itV2->second <= prev(endV1)->second)
    if !(a[0] <= b[b.len() - 1] || b[0] <= a[a.len() - 1]) {
        return;
    }

    let (mut i1, mut i2) = (0usize, 0usize);
    while i1 < a.len() && i2 < b.len() {
        if a[i1] < b[i2] {
            i1 += 1;
        } else if b[i2] < a[i1] {
            i2 += 1;
        } else {
            // if (!mipsolver.submip) printf("node %d can be pruned\n", itV2->second);
            let pruned_node = b[i2];
            i1 += 1;
            i2 += 1;
            let w = ms.mipdata_.nodequeue.prune_node(pruned_node);
            ms.mipdata_.pruned_treeweight += w;
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests

    use super::*;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;

    /// A clique table over `ncols` columns with a fixed pseudo-random set of cliques
    /// (sizes 2..=5, both literal polarities), so that some pairs share a clique and some
    /// do not, and some cliques have size two (`sizeTwoCliques`).
    fn build_table(ncols: HighsInt) -> HighsCliqueTable {
        let mut ct = HighsCliqueTable::new(ncols);
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as HighsInt
        };
        for _ in 0..40 {
            let len = 2 + next() % 4;
            let mut cols: Vec<HighsInt> = Vec::new();
            while (cols.len() as HighsInt) < len {
                let col = next() % ncols;
                if !cols.contains(&col) {
                    cols.push(col);
                }
            }
            let clq: Vec<HighsCliqueTable_CliqueVar> = cols
                .iter()
                .map(|&col| HighsCliqueTable_CliqueVar::new(col, next() % 2))
                .collect();
            ct.do_add_clique(&clq[..], false, kHighsIInf);
        }
        ct
    }

    /// Every literal of every column: `N` = 2 * ncols, well above the grain size 10.
    fn all_literals(ncols: HighsInt) -> Vec<HighsCliqueTable_CliqueVar> {
        let mut q = Vec::new();
        for col in 0..ncols {
            q.push(HighsCliqueTable_CliqueVar::new(col, 0));
            q.push(HighsCliqueTable_CliqueVar::new(col, 1));
        }
        q
    }

    /// Brute force over HighsCliqueTable.cpp 489-491: the indices `i` with
    /// `haveCommonClique(v, q[i])`, and the summed query count.
    fn brute_force(
        ct: &HighsCliqueTable,
        v: HighsCliqueTable_CliqueVar,
        q: &[HighsCliqueTable_CliqueVar],
    ) -> (Vec<HighsInt>, i64) {
        let mut expected = Vec::new();
        let mut nq: i64 = 0;
        for (i, qi) in q.iter().enumerate() {
            if ct.have_common_clique_2(&mut nq, v, *qi) {
                expected.push(i as HighsInt);
            }
        }
        (expected, nq)
    }

    #[test]
    fn w5_mip_clique_query_neighbourhood_parallel() {
        // `queryNeighbourhood` (HighsCliqueTable.cpp 480-520) runs `for_each` (grain 10) into a
        // `HighsCombinable` and `combine_each` + pdqsort give the serial neighbourhood and
        // `numNeighbourhoodQueries` at threads 1, 2, 4 (and 8), each run twice.
        let ncols: HighsInt = 60;
        let ct = build_table(ncols);
        let q = all_literals(ncols);
        let num_vars = q.len();
        assert!(
            ct.size_two_cliques.size() > 0,
            "fixture has size-two cliques"
        );

        // Serial branch (minEntriesForParallelism = kHighsIInf, the threads = 1 setting).
        let mut serial: Vec<(Vec<HighsInt>, i64)> = Vec::new();
        let mut any_nonempty = false;
        for k in 0..num_vars {
            let v = q[k];
            let mut inds: Vec<HighsInt> = vec![-1, -2];
            let mut nq: i64 = 7;
            ct.query_neighbourhood(&mut inds, &mut nq, v, &q[..]);
            if ct.num_cliques_2(v) == 0 {
                // C++ line 485: early return after the clear, numQueries untouched.
                assert!(inds.is_empty());
                assert_eq!(nq, 7);
            } else {
                let (expected, expected_nq) = brute_force(&ct, v, &q[..]);
                assert_eq!(inds, expected, "serial neighbourhood of literal {k}");
                assert_eq!(nq, 7 + expected_nq, "serial numQueries of literal {k}");
                any_nonempty |= !inds.is_empty();
            }
            serial.push((inds, nq));
        }
        assert!(any_nonempty, "fixture has non-empty neighbourhoods");

        // allowParallel = false takes the serial branch even with minEntriesForParallelism = 0
        // (C++ line 487).
        {
            let mut ct_np = ct.clone();
            ct_np.set_min_entries_for_parallelism(0);
            ct_np.set_allow_parallel(false);
            for k in 0..num_vars {
                let mut inds: Vec<HighsInt> = Vec::new();
                let mut nq: i64 = 7;
                ct_np.query_neighbourhood(&mut inds, &mut nq, q[k], &q[..]);
                assert_eq!((inds, nq), serial[k], "allowParallel=false, literal {k}");
            }
        }

        // Parallel branch: minEntriesForParallelism = 0, so `numEntries - 2 * |sizeTwo| < 0` is
        // false (C++ line 488). Each thread count runs on its own thread (own executor).
        for threads in [1, 2, 4, 8] {
            for _rep in 0..2 {
                let mut ct_par = ct.clone();
                ct_par.set_min_entries_for_parallelism(0);
                let q_par = q.clone();
                let serial_par = serial.clone();
                std::thread::spawn(move || {
                    initialize_scheduler(threads);
                    for k in 0..q_par.len() {
                        let mut inds: Vec<HighsInt> = vec![-1];
                        let mut nq: i64 = 7;
                        ct_par.query_neighbourhood(&mut inds, &mut nq, q_par[k], &q_par[..]);
                        assert_eq!(
                            (inds, nq),
                            serial_par[k],
                            "parallel branch at threads = {threads}, literal {k}"
                        );
                    }
                    HighsTaskExecutor::shutdown(true);
                })
                .join()
                .unwrap();
            }
        }
    }

    /// A clique table as `new_data` builds it with `h16_clique_links` set to `on`.
    fn h16_clique_links_table(ncols: HighsInt, on: bool) -> HighsCliqueTable {
        let mut opts = HighsOptions::default();
        opts.h16_clique_links = on;
        let mut ct = HighsCliqueTable::new(ncols);
        ct.set_h16_flags(&opts);
        ct
    }

    #[test]
    fn h16_clique_links_on_duplicate_size_two_cliques() {
        // check/TestMipSolver.cpp "MIP-clique-table-duplicate-size-two-cliques" (1.16.0 diff
        // lines 4-48): several size-two cliques on one edge keep the edge until the last one
        // is removed.
        let edge = [
            HighsCliqueTable_CliqueVar::new(0, 1),
            HighsCliqueTable_CliqueVar::new(1, 1),
        ];
        let other_edge = [
            HighsCliqueTable_CliqueVar::new(0, 1),
            HighsCliqueTable_CliqueVar::new(2, 1),
        ];

        // remove the clique added last
        {
            let mut ct = h16_clique_links_table(3, true);
            ct.do_add_clique(&edge, false, kHighsIInf);
            ct.do_add_clique(&edge, false, kHighsIInf);
            assert_eq!(ct.num_cliques(), 2);
            ct.remove_clique(1, true);
            assert!(ct.have_common_clique(edge[0], edge[1]));
            ct.remove_clique(0, true);
            assert!(!ct.have_common_clique(edge[0], edge[1]));
        }

        // remove the clique added first
        {
            let mut ct = h16_clique_links_table(3, true);
            ct.do_add_clique(&edge, false, kHighsIInf);
            ct.do_add_clique(&other_edge, false, kHighsIInf);
            ct.do_add_clique(&edge, false, kHighsIInf);
            ct.remove_clique(0, true);
            assert!(ct.have_common_clique(edge[0], edge[1]));
            assert!(ct.have_common_clique(other_edge[0], other_edge[1]));
            ct.remove_clique(2, true);
            assert!(!ct.have_common_clique(edge[0], edge[1]));
            assert!(ct.have_common_clique(other_edge[0], other_edge[1]));
        }
    }

    #[test]
    fn h16_clique_links_on_substitute_and_remove() {
        // replaceLiteral (1.16.0 HighsCliqueTable.cpp lines 641-676): both size-two cliques on
        // the edge move to (0,1)-(3,1) and stay linked on that edge until both are removed.
        let mut ct = h16_clique_links_table(4, true);
        let edge = [
            HighsCliqueTable_CliqueVar::new(0, 1),
            HighsCliqueTable_CliqueVar::new(1, 1),
        ];
        ct.do_add_clique(&edge, false, kHighsIInf);
        ct.do_add_clique(&edge, false, kHighsIInf);
        ct.replace_substituted_lists(
            HighsCliqueTable_CliqueVar::new(1, 1),
            HighsCliqueTable_CliqueVar::new(3, 1),
        );
        assert!(!ct.have_common_clique(edge[0], edge[1]));
        assert!(ct.have_common_clique(
            HighsCliqueTable_CliqueVar::new(0, 1),
            HighsCliqueTable_CliqueVar::new(3, 1)
        ));
        ct.remove_clique(0, true);
        assert!(ct.have_common_clique(
            HighsCliqueTable_CliqueVar::new(0, 1),
            HighsCliqueTable_CliqueVar::new(3, 1)
        ));
        ct.remove_clique(1, true);
        assert!(!ct.have_common_clique(
            HighsCliqueTable_CliqueVar::new(0, 1),
            HighsCliqueTable_CliqueVar::new(3, 1)
        ));
    }

    #[test]
    fn h16_clique_links_off_matches_1_15_1_table() {
        // The same sequences with the option off reproduce the 1.15.1 table: removeClique and
        // replaceLiteral erase the edge key outright (1.15.1 HighsCliqueTable.cpp lines 823 and
        // 631), so the edge is gone while another size-two clique still holds it.
        let edge = [
            HighsCliqueTable_CliqueVar::new(0, 1),
            HighsCliqueTable_CliqueVar::new(1, 1),
        ];
        let other_edge = [
            HighsCliqueTable_CliqueVar::new(0, 1),
            HighsCliqueTable_CliqueVar::new(2, 1),
        ];

        {
            let mut ct = h16_clique_links_table(3, false);
            ct.do_add_clique(&edge, false, kHighsIInf);
            ct.do_add_clique(&edge, false, kHighsIInf);
            assert_eq!(ct.num_cliques(), 2);
            assert!(ct.size_two_links.is_empty());
            ct.remove_clique(1, true);
            assert!(!ct.have_common_clique(edge[0], edge[1]));
        }

        {
            let mut ct = h16_clique_links_table(3, false);
            ct.do_add_clique(&edge, false, kHighsIInf);
            ct.do_add_clique(&other_edge, false, kHighsIInf);
            ct.do_add_clique(&edge, false, kHighsIInf);
            ct.remove_clique(0, true);
            assert!(!ct.have_common_clique(edge[0], edge[1]));
            assert!(ct.have_common_clique(other_edge[0], other_edge[1]));
            ct.remove_clique(2, true);
            assert!(!ct.have_common_clique(edge[0], edge[1]));
            assert!(ct.have_common_clique(other_edge[0], other_edge[1]));
        }

        {
            let mut ct = h16_clique_links_table(4, false);
            ct.do_add_clique(&edge, false, kHighsIInf);
            ct.do_add_clique(&edge, false, kHighsIInf);
            ct.replace_substituted_lists(
                HighsCliqueTable_CliqueVar::new(1, 1),
                HighsCliqueTable_CliqueVar::new(3, 1),
            );
            let moved = [
                HighsCliqueTable_CliqueVar::new(0, 1),
                HighsCliqueTable_CliqueVar::new(3, 1),
            ];
            assert!(ct.have_common_clique(moved[0], moved[1]));
            ct.remove_clique(1, true);
            assert!(!ct.have_common_clique(moved[0], moved[1]));
        }
    }

    /// HiGHS backport-AGG-B: `cliqueCover` (HighsCliqueTable.cpp 1082-1134) covers every literal, and the
    /// removal pass leaves no redundant clique: each kept clique has a member that occurs in no
    /// other kept clique.
    #[test]
    fn x1_agg_b_clique_cover_drops_redundant() {
        let mut ct = build_table(12);
        let mut clq_vars = all_literals(12);
        let mut cover: Vec<Vec<HighsCliqueTable_CliqueVar>> = Vec::new();
        ct.clique_cover(&mut clq_vars, &mut cover);
        assert!(!cover.is_empty());

        for lit in all_literals(12) {
            assert!(
                cover.iter().any(|clique| clique.contains(&lit)),
                "literal not covered"
            );
        }

        let mut count: HighsHashTable<HighsInt, HighsInt> = HighsHashTable::new();
        for clique in cover.iter() {
            for var in clique.iter() {
                *count.index_or_insert(&var.index()) += 1;
            }
        }
        for clique in cover.iter() {
            let has_unique = clique
                .iter()
                .any(|var| *count.index_or_insert(&var.index()) == 1);
            assert!(has_unique, "redundant clique kept in cover");
        }
    }
}
