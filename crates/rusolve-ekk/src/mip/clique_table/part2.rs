//! HiGHS source: mip/HighsCliqueTable.cpp (lines 1077-2243)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;
use crate::util::cpp_std::cpp_partition;
use crate::util::highs_pdqsort::pdqsort_by;
use crate::util::highs_utils::{cpp_max, cpp_min};

thread_local! {
    /// Rust-only scratch (perf): clique ids collected by `add_implications`, kept to reuse
    /// its allocation. A stack shared by the nested calls of one thread (the function
    /// re-enters itself through `change_bound`); empty between top-level calls. It is per
    /// thread, not a field of the clique table: parallel MIP workers share the table (it is
    /// only read here, as in C++) and each runs its own `add_implications` concurrently.
    /// Never borrowed across `change_bound`.
    static ADD_IMPLICATIONS_IDS: std::cell::RefCell<Vec<HighsInt>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

impl HighsCliqueTable {
    /// HiGHS 1.16 #3540 (h16_clique_links): option plumbing with no C++ counterpart. Reads the option
    /// once; the flag is fixed at table construction (see `add_size_two_clique`, 1.16.0
    /// HighsCliqueTable.cpp lines 34-51).
    pub fn set_h16_flags(&mut self, opts: &HighsOptions) {
        self.h16_clique_links = opts.h16_clique_links;
        // HiGHS 1.16 #3327 (h16_hash_common): option plumbing, read by find_common_clique_id
        self.h16_hash_common = opts.h16_hash_common;
        // HiGHS #3312 (mip_presolve_clique_table): option plumbing, the clique-table gate
        self.x1b_presolve_clique_table = opts.mip_presolve_clique_table;
    }

    /// HiGHS 1.16 #3540 (h16_clique_links): option plumbing with no C++ counterpart. Copies the flag
    /// into a table that `rebuild` or `build_from` creates from this one.
    pub fn copy_h16_flags_from(&mut self, other: &HighsCliqueTable) {
        self.h16_clique_links = other.h16_clique_links;
        // HiGHS 1.16 #3327 (h16_hash_common): option plumbing, copied with the other HiGHS 1.16 flags
        self.h16_hash_common = other.h16_hash_common;
        // HiGHS #3312 (mip_presolve_clique_table): option plumbing, copied with the HiGHS 1.16 flags
        self.x1b_presolve_clique_table = other.x1b_presolve_clique_table;
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::extractCliquesFromCut (lines 1078-1266)
    /// C++: `void extractCliquesFromCut(const HighsMipSolver & mipsolver, const HighsInt * inds, const double * vals, HighsInt len, double rhs)`
    /// PORT NOTE: signature changed to `ms: &mut` because the C++ `const HighsMipSolver&` mutates the
    /// global domain and the implications through `mipdata_` (M1 ruling).
    #[allow(clippy::needless_late_init)]
    pub fn extract_cliques_from_cut(
        &mut self,
        ms: &mut HighsMipSolver,
        inds: &[HighsInt],
        vals: &[f64],
        len: HighsInt,
        rhs: f64,
    ) {
        // HiGHS #3223 (mip_cdouble_proofs): HighsCliqueTable.cpp lines 1190-1280 at 7db1c8180
        let x1b_cdouble = ms.options_mip_.mip_cdouble_proofs;

        if self.is_full() {
            return;
        }

        let feastol = ms.mipdata_.feastol;

        let mut minact = HighsCDouble::from(0.0);
        let mut nbin: HighsInt = 0;
        for i in 0..len as usize {
            if ms.mipdata_.domains[0].is_binary(ms, inds[i]) {
                nbin += 1;
            }

            if vals[i] > 0.0 {
                if ms.mipdata_.domains[0].col_lower_[inds[i] as usize] == -kHighsInf {
                    return;
                }
                // HiGHS #3223 (mip_cdouble_proofs): HighsCliqueTable.cpp lines 1190-1191 at 7db1c8180
                if x1b_cdouble {
                    minact += vals[i]
                        * HighsCDouble::from(ms.mipdata_.domains[0].col_lower_[inds[i] as usize]);
                } else {
                    minact += vals[i] * ms.mipdata_.domains[0].col_lower_[inds[i] as usize];
                }
            } else {
                if ms.mipdata_.domains[0].col_upper_[inds[i] as usize] == kHighsInf {
                    return;
                }
                // HiGHS #3223 (mip_cdouble_proofs): HighsCliqueTable.cpp lines 1194-1195 at 7db1c8180
                if x1b_cdouble {
                    minact += vals[i]
                        * HighsCDouble::from(ms.mipdata_.domains[0].col_upper_[inds[i] as usize]);
                } else {
                    minact += vals[i] * ms.mipdata_.domains[0].col_upper_[inds[i] as usize];
                }
            }
        }

        if rhs - minact < 0.0 {
            minact = HighsCDouble::from(rhs);
        }

        // PORT NOTE (perf): the table is lent back to its slot once around the bound-change loop
        // (for `changeBound` -> `addImplications`) instead of once per `changeBound` (two swaps of
        // the whole `HighsCliqueTable` per integer column of every cut); the loop never reads `self`.
        let lent = ms.mipdata_.cliquetable.lend(self);
        let infeasible = 'bounds: {
            for i in 0..len as usize {
                if ms.is_col_continuous(inds[i]) {
                    continue;
                }

                // HiGHS #3223 (mip_cdouble_proofs): HighsCliqueTable.cpp lines 1206-1221 at 7db1c8180
                if x1b_cdouble {
                    let implied_bound: HighsCDouble = (rhs - minact) / vals[i];
                    if vals[i] > 0.0 {
                        let domfeas = ms.mipdata_.domains[0].feastol(ms);
                        let bound_tol = cpp_max(domfeas, ms.mipdata_.epsilon / vals[i].abs());
                        let bound_val = f64::from(
                            implied_bound
                                + ms.mipdata_.domains[0].col_lower_[inds[i] as usize]
                                + bound_tol,
                        )
                        .floor();
                        // PORT NOTE: C++ `changeBound(kUpper, col, val, reason)` is the struct form here.
                        detach!(ms.mipdata_.domains[0], |d| d.change_bound(
                            ms,
                            HighsDomainChange {
                                boundval: bound_val,
                                column: inds[i],
                                boundtype: HighsBoundType::kUpper,
                            },
                            HighsDomain_Reason::unspecified()
                        ));
                        if ms.mipdata_.domains[0].infeasible() {
                            break 'bounds true;
                        }
                    } else {
                        let domfeas = ms.mipdata_.domains[0].feastol(ms);
                        let bound_tol = cpp_max(domfeas, ms.mipdata_.epsilon / vals[i].abs());
                        let bound_val = f64::from(
                            implied_bound + ms.mipdata_.domains[0].col_upper_[inds[i] as usize]
                                - bound_tol,
                        )
                        .ceil();
                        detach!(ms.mipdata_.domains[0], |d| d.change_bound(
                            ms,
                            HighsDomainChange {
                                boundval: bound_val,
                                column: inds[i],
                                boundtype: HighsBoundType::kLower,
                            },
                            HighsDomain_Reason::unspecified()
                        ));
                        if ms.mipdata_.domains[0].infeasible() {
                            break 'bounds true;
                        }
                    }
                } else {
                    let mut bound_val = f64::from((rhs - minact) / vals[i]);
                    if vals[i] > 0.0 {
                        let domfeas = ms.mipdata_.domains[0].feastol(ms);
                        bound_val = (bound_val
                            + ms.mipdata_.domains[0].col_lower_[inds[i] as usize]
                            + domfeas)
                            .floor();
                        // PORT NOTE: C++ `changeBound(kUpper, col, val, reason)` is the struct form here.
                        detach!(ms.mipdata_.domains[0], |d| d.change_bound(
                            ms,
                            HighsDomainChange {
                                boundval: bound_val,
                                column: inds[i],
                                boundtype: HighsBoundType::kUpper,
                            },
                            HighsDomain_Reason::unspecified()
                        ));
                        if ms.mipdata_.domains[0].infeasible() {
                            break 'bounds true;
                        }
                    } else {
                        let domfeas = ms.mipdata_.domains[0].feastol(ms);
                        bound_val = (bound_val
                            + ms.mipdata_.domains[0].col_upper_[inds[i] as usize]
                            - domfeas)
                            .ceil();
                        detach!(ms.mipdata_.domains[0], |d| d.change_bound(
                            ms,
                            HighsDomainChange {
                                boundval: bound_val,
                                column: inds[i],
                                boundtype: HighsBoundType::kLower,
                            },
                            HighsDomain_Reason::unspecified()
                        ));
                        if ms.mipdata_.domains[0].infeasible() {
                            break 'bounds true;
                        }
                    }
                }
            }
            false
        };
        ms.mipdata_.cliquetable.reclaim(self, lent);
        if infeasible {
            return;
        }

        if nbin <= 1 {
            return;
        }

        let mut perm: Vec<HighsInt> = (0..len).collect();

        let globaldom = &ms.mipdata_.domains[0];
        let binaryend = cpp_partition(&mut perm, |pos: &HighsInt| {
            globaldom.is_binary(ms, inds[*pos as usize])
        });

        nbin = binaryend as HighsInt;

        // if not all variables are binary, we extract variable upper and lower bounds
        // constraints on the non-binary variable for each binary variable in the
        // constraint:
        if nbin < len {
            for i in 0..nbin as usize {
                let bincol = inds[perm[i] as usize];
                let mut implied_activity = rhs - minact - vals[perm[i] as usize].abs();
                for j in nbin as usize..len as usize {
                    let col = inds[perm[j] as usize];
                    if ms.mipdata_.domains[0].is_fixed(col) {
                        continue;
                    }

                    if vals[perm[j] as usize] > 0.0 {
                        // HiGHS #3223 (mip_cdouble_proofs): HighsCliqueTable.cpp lines 1248-1252 at 7db1c8180
                        let mut implcolub = if x1b_cdouble {
                            f64::from(
                                (implied_activity
                                    + vals[perm[j] as usize]
                                        * HighsCDouble::from(
                                            ms.mipdata_.domains[0].col_lower_[col as usize],
                                        ))
                                    / vals[perm[j] as usize],
                            )
                        } else {
                            // C++ 1.15.1 `static_cast<double>(impliedActivity + a * bound) / a`:
                            // the cast binds to the sum, so the division is in double.
                            f64::from(
                                implied_activity
                                    + vals[perm[j] as usize]
                                        * ms.mipdata_.domains[0].col_lower_[col as usize],
                            ) / vals[perm[j] as usize]
                        };
                        if ms.is_col_integral(col) {
                            implcolub = (implcolub + ms.mipdata_.feastol).floor();
                        }

                        if implcolub < ms.mipdata_.domains[0].col_upper_[col as usize] - feastol {
                            let coef;
                            let constant;
                            if vals[perm[i] as usize] < 0.0 {
                                coef = ms.mipdata_.domains[0].col_upper_[col as usize] - implcolub;
                                constant = implcolub;
                            } else {
                                // make sure that upper bound is not infinite to avoid adding VUB
                                // with coefficient '-kHighsInf' and constant 'kHighsInf'
                                if ms.mipdata_.domains[0].col_upper_[col as usize] == kHighsInf {
                                    continue;
                                }
                                coef = implcolub - ms.mipdata_.domains[0].col_upper_[col as usize];
                                constant = ms.mipdata_.domains[0].col_upper_[col as usize];
                            }
                            detach!(ms.mipdata_.implications, |im| im
                                .add_vub(ms, col, bincol, coef, constant));
                        }
                    } else {
                        // HiGHS #3223 (mip_cdouble_proofs): HighsCliqueTable.cpp lines 1275-1278 at 7db1c8180
                        let mut implcollb = if x1b_cdouble {
                            f64::from(
                                (implied_activity
                                    + vals[perm[j] as usize]
                                        * HighsCDouble::from(
                                            ms.mipdata_.domains[0].col_upper_[col as usize],
                                        ))
                                    / vals[perm[j] as usize],
                            )
                        } else {
                            // C++ 1.15.1 `static_cast<double>(impliedActivity + a * bound) / a`:
                            // the cast binds to the sum, so the division is in double.
                            f64::from(
                                implied_activity
                                    + vals[perm[j] as usize]
                                        * ms.mipdata_.domains[0].col_upper_[col as usize],
                            ) / vals[perm[j] as usize]
                        };
                        if ms.is_col_integral(col) {
                            implcollb = (implcollb - ms.mipdata_.feastol).ceil();
                        }

                        if implcollb > ms.mipdata_.domains[0].col_lower_[col as usize] + feastol {
                            let coef;
                            let constant;
                            if vals[perm[i] as usize] < 0.0 {
                                coef = ms.mipdata_.domains[0].col_lower_[col as usize] - implcollb;
                                constant = implcollb;
                            } else {
                                // make sure that lower bound is not infinite to avoid adding VLB
                                // with coefficient 'kHighsInf' and constant '-kHighsInf'
                                if ms.mipdata_.domains[0].col_lower_[col as usize] == -kHighsInf {
                                    continue;
                                }
                                coef = implcollb - ms.mipdata_.domains[0].col_lower_[col as usize];
                                constant = ms.mipdata_.domains[0].col_lower_[col as usize];
                            }
                            detach!(ms.mipdata_.implications, |im| im
                                .add_vlb(ms, col, bincol, coef, constant));
                        }
                    }
                }
                let _ = &mut implied_activity;
            }
        }

        // only one binary means we do have no cliques
        if nbin <= 1 {
            return;
        }

        let mut clique: Vec<HighsCliqueTable_CliqueVar> = Vec::with_capacity(nbin as usize);

        pdqsort_by(
            &mut perm[..nbin as usize],
            |p1: &HighsInt, p2: &HighsInt| {
                // std::make_pair(std::abs(vals[p1]), p1) > std::make_pair(std::abs(vals[p2]), p2)
                let a1 = vals[*p1 as usize].abs();
                let a2 = vals[*p2 as usize].abs();
                a1 > a2 || (a1 == a2 && *p1 > *p2)
            },
        );
        // check if any cliques exists
        if vals[perm[0] as usize].abs() + vals[perm[1] as usize].abs()
            <= f64::from(rhs - minact + feastol)
        {
            return;
        }

        let nnz = ms.mipdata_.domains[0].num_model_nonzeros(ms);
        let max_new_entries = cpp_min(
            ms.mipdata_.num_clique_entries_after_presolve + 100000 + 4 * nnz,
            self.num_entries + 10 * nbin,
        );

        // PORT NOTE: `for (k = nbin - 1; k != 0 && numEntries < maxNewEntries; --k)` written as a
        // descending range with the same loop-condition check at the top of each iteration.
        for k in (1..nbin as usize).rev() {
            if self.num_entries >= max_new_entries {
                break;
            }

            let mincliqueval = f64::from(rhs - minact - vals[perm[k] as usize].abs() + feastol);
            let cliqueend = perm[..k].partition_point(|p| vals[*p as usize].abs() > mincliqueval);

            // no clique for this variable
            if cliqueend == 0 {
                continue;
            }

            clique.clear();

            for j in 0..cliqueend {
                let pos = perm[j] as usize;
                if vals[pos] < 0.0 {
                    clique.push(HighsCliqueTable_CliqueVar::new(inds[pos], 0));
                } else {
                    clique.push(HighsCliqueTable_CliqueVar::new(inds[pos], 1));
                }
            }

            let kp = perm[k] as usize;
            if vals[kp] < 0.0 {
                clique.push(HighsCliqueTable_CliqueVar::new(inds[kp], 0));
            } else {
                clique.push(HighsCliqueTable_CliqueVar::new(inds[kp], 1));
            }

            // printf("extracted this clique:\n");
            // printClique(clique);
            if clique.len() >= 2 {
                // PORT NOTE: `addClique(mipsolver, clique.data(), size)` as slice (M1 ruling).
                self.add_clique(ms, &mut clique, false, kHighsIInf);
                if ms.mipdata_.domains[0].infeasible() || self.num_entries >= max_new_entries {
                    return;
                }
            }

            // further cliques are just subsets of this clique
            if cliqueend == k {
                return;
            }
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::extractCliques (lines 1268-1378)
    /// HiGHS latest (#3286, 843e9fdb6): HighsCliqueTable.cpp::HighsCliqueTable::extractCliques (lines 1295-1360)
    /// C++: `void extractCliques(HighsMipSolver & mipsolver, bool transformRows = true)`
    /// C++ default arguments: `transform_rows = true`
    #[allow(clippy::needless_return)] // C++ `if (globaldom.infeasible()) return;` kept
    pub fn extract_cliques_2(&mut self, ms: &mut HighsMipSolver, transform_rows: bool) {
        let mut inds: Vec<HighsInt> = Vec::new();
        let mut vals: Vec<f64> = Vec::new();
        let mut perm: Vec<HighsInt> = Vec::new();
        let mut complementation: Vec<i8> = Vec::new();
        let mut clique: Vec<HighsCliqueTable_CliqueVar> = Vec::new();
        let mut entries: HighsHashTable<HighsInt, f64> = HighsHashTable::new();
        let mut offset: f64;

        // HiGHS #3286 (mip_clique_normalise): HighsCliqueTable.cpp lines 1295-1360 at 843e9fdb6
        let x1_clique_normalise = ms.options_mip_.mip_clique_normalise;

        // HiGHS #3250 (mip_presolve_row_append): HighsCliqueTable.cpp line 1414 (1.16.0)
        // loops over getNonCutRows(), so rows appended by FM presolve are scanned too.
        let x1b_row_append = ms.options_mip_.mip_presolve_row_append;
        let rows: Vec<HighsInt> = if x1b_row_append {
            ms.mipdata_.post_solve_stack.get_non_cut_rows()
        } else {
            (0..ms.num_row()).collect()
        };
        for i in rows {
            let iu = i as usize;
            let start = ms.mipdata_.a_rstart_[iu];
            let end = ms.mipdata_.a_rstart_[iu + 1];

            if !x1b_row_append
                && ms.mipdata_.post_solve_stack.get_orig_row_index(i) >= ms.orig_model_.num_row_
            {
                break;
            }

            // catch set packing and partitioning constraints that already have the form
            // of a clique without transformations and add those cliques with the rows
            // being recorded
            if x1_clique_normalise {
                // catch set packing and partitioning constraints that already have the
                // form of a clique without transformations and add those cliques with
                // the rows being recorded. only <= and = rows are checked because
                // normaliseCliqueRows has already flipped >= rows to <= form.
                if ms.row_upper_2(i) < kHighsInf {
                    let mut issetppc = true;
                    let mut fixed_term = HighsCDouble::from(0.0);
                    let mut num_comp: HighsInt = 0;
                    for j in start..end {
                        let col = ms.mipdata_.a_rindex_[j as usize];
                        let val = ms.mipdata_.a_rvalue_[j as usize];

                        // handle fixed non-binary variables
                        if !detach!(ms.mipdata_.domains[0], |d| d.is_binary(ms, col))
                            && ms.mipdata_.domains[0].is_fixed(col)
                        {
                            fixed_term += val
                                * HighsCDouble::from(
                                    ms.mipdata_.domains[0].col_upper_[col as usize],
                                );
                            continue;
                        }

                        // check if we have a set partitioning / packing row
                        issetppc = detach!(ms.mipdata_.domains[0], |d| d.is_binary(ms, col))
                            && val.abs() == 1.0;
                        if !issetppc {
                            break;
                        }

                        // count number of complemented binaries (with coefficient -1)
                        if val < 0.0 {
                            num_comp += 1;
                        }
                    }

                    if issetppc {
                        // subtract fixed term and round right-hand side
                        let rhs = f64::from(
                            (ms.row_upper_2(i) - fixed_term + ms.mipdata_.feastol).floor(),
                        );

                        if rhs == 1.0 - num_comp as f64 {
                            // subtract fixed term and round left-hand side if finite
                            let mut lhs = ms.row_lower_2(i);
                            if lhs > -kHighsInf {
                                lhs = f64::from((lhs - fixed_term - ms.mipdata_.feastol).ceil());
                            }

                            clique.clear();

                            for j in start..end {
                                let col = ms.mipdata_.a_rindex_[j as usize];
                                let val = ms.mipdata_.a_rvalue_[j as usize];
                                let dir: HighsInt = if val > 0.0 { 1 } else { 0 };

                                // skip non-binary variables (fixed, see previous loop) and binaries
                                // that are fixed to "inactive" values
                                if !detach!(ms.mipdata_.domains[0], |d| d.is_binary(ms, col))
                                    || ms.mipdata_.domains[0].is_fixed_to_val(col, 1 - dir)
                                {
                                    continue;
                                }

                                // add to clique
                                clique.push(HighsCliqueTable_CliqueVar::new(col, dir));
                            }

                            // add clique to clique table
                            if clique.len() >= 2 {
                                self.add_clique(ms, &mut clique, rhs == lhs, i);
                                if ms.mipdata_.domains[0].infeasible() {
                                    return;
                                }
                            }
                            continue;
                        }
                    }
                }
            } else {
                // catch set packing and partitioning constraints that already have the form
                // of a clique without transformations and add those cliques with the rows
                // being recorded
                if ms.row_upper_2(i) == 1.0 {
                    let mut issetppc = true;

                    clique.clear();

                    for j in start..end {
                        let col = ms.mipdata_.a_rindex_[j as usize];
                        if ms.mipdata_.domains[0].col_upper_[col as usize] == 0.0
                            && ms.mipdata_.domains[0].col_lower_[col as usize] == 0.0
                        {
                            continue;
                        }

                        issetppc = detach!(ms.mipdata_.domains[0], |d| d.is_binary(ms, col))
                            && ms.mipdata_.a_rvalue_[j as usize] == 1.0;
                        if !issetppc {
                            break;
                        }

                        clique.push(HighsCliqueTable_CliqueVar::new(col, 1));
                    }

                    if issetppc {
                        let equality = ms.row_lower_2(i) == 1.0;
                        self.add_clique(ms, &mut clique, equality, i);
                        if ms.mipdata_.domains[0].infeasible() {
                            return;
                        }
                        continue;
                    }
                }
            }
            if !transform_rows || self.is_full() {
                continue;
            }

            offset = 0.0;
            for j in start..end {
                let mut col = ms.mipdata_.a_rindex_[j as usize];
                let mut val = ms.mipdata_.a_rvalue_[j as usize];

                self.resolve_substitution_2(&mut col, &mut val, &mut offset);
                *entries.index_or_insert(&col) += val;
            }

            // PORT NOTE: the row bounds are read before the `checkRow` closure is built, because
            // the closure borrows `ms` mutably. `checkRow` does not change the model rows, so
            // the values equal the ones C++ reads at each call.
            let row_upper = ms.row_upper_2(i);
            let row_lower = ms.row_lower_2(i);
            let mut check_row = |rhs0: f64, direction: HighsInt| {
                // if (direction * rhs == kHighsInf) return;
                if direction as f64 * rhs0 == kHighsInf {
                    return;
                }
                let mut rhs = direction as f64 * (rhs0 - offset);
                inds.clear();
                vals.clear();
                complementation.clear();
                let mut freevar = false;
                let mut nbin: HighsInt = 0;

                for entry in entries.iter() {
                    let col = *entry.key();
                    let val = direction as f64 * *entry.value();

                    if val.abs() < ms.mipdata_.epsilon {
                        continue;
                    }

                    if detach!(ms.mipdata_.domains[0], |d| d.is_binary(ms, col)) {
                        nbin += 1;
                    }

                    if val < 0.0 {
                        freevar = ms.mipdata_.domains[0].col_upper_[col as usize] == kHighsInf;
                        if freevar {
                            break;
                        }

                        vals.push(-val);
                        inds.push(col);
                        complementation.push(-1);
                        rhs -= val * ms.mipdata_.domains[0].col_upper_[col as usize];
                    } else {
                        freevar = ms.mipdata_.domains[0].col_lower_[col as usize] == -kHighsInf;
                        if freevar {
                            break;
                        }

                        vals.push(val);
                        inds.push(col);
                        complementation.push(1);
                        rhs -= val * ms.mipdata_.domains[0].col_lower_[col as usize];
                    }
                }

                if !freevar && nbin != 0 {
                    let feastol = ms.mipdata_.feastol;
                    self.extract_cliques(
                        ms,
                        &mut inds,
                        &mut vals,
                        &mut complementation,
                        rhs,
                        nbin,
                        &mut perm,
                        &mut clique,
                        feastol,
                    );
                    if ms.mipdata_.domains[0].infeasible() {
                        return;
                    }
                }
            };

            check_row(row_upper, 1);
            check_row(row_lower, -1);

            entries.clear();
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::extractObjCliques (lines 1380-1461)
    /// C++: `void extractObjCliques(HighsMipSolver & mipsolver)`
    pub fn extract_obj_cliques(&mut self, ms: &mut HighsMipSolver) {
        let mut nbin = ms
            .mipdata_
            .objective_function
            .get_num_binaries_in_objective();
        if nbin <= 1 {
            return;
        }
        if ms.mipdata_.domains[0].get_objective_lower_bound() == -kHighsInf {
            return;
        }

        // PORT NOTE: `getCutoffConstraint` hands out slices borrowed from the domain. They are
        // copied, so the domain can be read while they are alive. The cutoff row is not changed
        // inside this function.
        let (cvals, cinds, clen, crhs) = detach!(ms.mipdata_.domains[0], |d| {
            let mut vals: &[f64] = &[];
            let mut inds: &[HighsInt] = &[];
            let mut len: HighsInt = 0;
            let mut rhs: f64 = 0.0;
            d.get_cutoff_constraint(ms, &mut vals, &mut inds, &mut len, &mut rhs);
            (vals.to_vec(), inds.to_vec(), len, rhs)
        });

        let mut perm: Vec<HighsInt> = (0..nbin).collect();

        let binaryend = cpp_partition(&mut perm, |pos: &HighsInt| {
            cvals[*pos as usize] != 0.0 && !ms.mipdata_.domains[0].is_fixed(cinds[*pos as usize])
        });

        nbin = binaryend as HighsInt;

        // only one binary means we do have no cliques
        if nbin <= 1 {
            return;
        }

        let mut clique: Vec<HighsCliqueTable_CliqueVar> = Vec::with_capacity(nbin as usize);

        pdqsort_by(
            &mut perm[..nbin as usize],
            |p1: &HighsInt, p2: &HighsInt| {
                let a1 = cvals[*p1 as usize].abs();
                let a2 = cvals[*p2 as usize].abs();
                a1 > a2 || (a1 == a2 && *p1 > *p2)
            },
        );

        // check if any cliques exists
        let mut minact = HighsCDouble::from(0.0);
        let mut ninf: HighsInt = 0;
        ms.mipdata_.domains[0].compute_min_activity(
            0,
            clen,
            &cinds,
            &cvals,
            &mut ninf,
            &mut minact,
        );
        let feastol = ms.mipdata_.feastol;
        if cvals[perm[0] as usize].abs() + cvals[perm[1] as usize].abs()
            <= f64::from(crhs - minact + feastol)
        {
            return;
        }

        for k in (1..nbin as usize).rev() {
            let mincliqueval = f64::from(crhs - minact - cvals[perm[k] as usize].abs() + feastol);
            let cliqueend = perm[..k].partition_point(|p| cvals[*p as usize].abs() > mincliqueval);

            // no clique for this variable
            if cliqueend == 0 {
                continue;
            }

            clique.clear();

            for j in 0..cliqueend {
                let pos = perm[j] as usize;
                if cvals[pos] < 0.0 {
                    clique.push(HighsCliqueTable_CliqueVar::new(cinds[pos], 0));
                } else {
                    clique.push(HighsCliqueTable_CliqueVar::new(cinds[pos], 1));
                }
            }

            let kp = perm[k] as usize;
            if cvals[kp] < 0.0 {
                clique.push(HighsCliqueTable_CliqueVar::new(cinds[kp], 0));
            } else {
                clique.push(HighsCliqueTable_CliqueVar::new(cinds[kp], 1));
            }

            if clique.len() >= 2 {
                self.add_clique(ms, &mut clique, false, kHighsIInf);
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

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::processInfeasibleVertices (lines 1463-1567)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::processInfeasibleVertices (lines 1591-1700)
    /// C++: `void processInfeasibleVertices(HighsDomain & globaldom)`
    /// PORT NOTE: signature changed to take `ms: &mut` (`fixCol` needs it, M1 ruling).
    pub fn process_infeasible_vertices(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldom: &mut HighsDomain,
    ) {
        // HiGHS #3220 (mip_fix_last_active): HighsCliqueTable.cpp lines 1533-1534, 1549 at b7d4e98c5
        let x1b_fix_last = ms.options_mip_.mip_fix_last_active;
        // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 1591-1700 at 7b890c255
        let clqp = self.x1b_presolve_clique_table;

        while !self.infeasvertexstack.is_empty() && !globaldom.infeasible() {
            let top = *self.infeasvertexstack.last().unwrap();
            self.infeasvertexstack.pop();
            let mut v = top.complement();

            self.resolve_substitution(&mut v);
            let wasfixed = globaldom.is_fixed(v.col() as HighsInt);
            lend!(
                ms.mipdata_.cliquetable,
                self,
                globaldom.fix_col(
                    ms,
                    v.col() as HighsInt,
                    v.val() as f64,
                    HighsDomain_Reason::unspecified(),
                )
            );
            if globaldom.infeasible() {
                return;
            }
            if !wasfixed {
                self.nfixings += 1;
            }
            if self.col_deleted[v.col() as usize] != 0 {
                continue;
            }
            self.col_deleted[v.col() as usize] = 1;

            let mut v_hash_lists = std::mem::take(&mut self.inverted_hash_list[v.index() as usize]);
            let mut v_hash_lists_size_two =
                std::mem::take(&mut self.inverted_hash_list_size_two[v.index() as usize]);

            let mut fix_clique = |cliqueid: HighsInt| -> bool {
                let start = self.cliques[cliqueid as usize].start;
                let end = self.cliques[cliqueid as usize].end;

                for i in start..end {
                    let ce = self.cliqueentries[i as usize];
                    if clqp {
                        // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 1614-1615, 1630-1631 at 7b890c255
                        if ce.col() == v.col() || self.col_deleted[ce.col() as usize] != 0 {
                            continue;
                        }
                    } else if ce.col() == v.col() {
                        continue;
                    }
                    if !self.fix_col(ms, globaldom, ce, false) {
                        return true;
                    }
                }

                self.remove_clique(cliqueid, true);
                false
            };

            let infeas =
                v_hash_lists.for_each_until(&mut |cid: &HighsInt, _: &HighsInt| fix_clique(*cid));
            if infeas {
                return;
            }

            let infeas = v_hash_lists_size_two
                .for_each_until(&mut |cid: &HighsInt, _: &()| fix_clique(*cid));
            if infeas {
                return;
            }

            v_hash_lists =
                std::mem::take(&mut self.inverted_hash_list[v.complement().index() as usize]);
            v_hash_lists_size_two = std::mem::take(
                &mut self.inverted_hash_list_size_two[v.complement().index() as usize],
            );

            if clqp {
                // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 1645-1697 at 7b890c255
                v_hash_lists_size_two.for_each(&mut |cliqueid: &HighsInt, _: &()| {
                    self.cliques[*cliqueid as usize].num_zero_fixed += 1;
                    self.fix_last_active_and_remove(ms, globaldom, *cliqueid);
                });

                debug_assert!(self.cliquehitinds.is_empty());
                let mut clq: Vec<HighsCliqueTable_CliqueVar> = Vec::new();
                v_hash_lists.for_each(&mut |cliqueid: &HighsInt, _: &HighsInt| {
                    let cid = *cliqueid as usize;
                    // assert(cliqueentries[entry.value()].val == 1 - v.val);

                    self.cliques[cid].num_zero_fixed += 1;
                    if self.cliques[cid].num_active() <= 1 {
                        self.fix_last_active_and_remove(ms, globaldom, *cliqueid);
                    } else if !self.in_presolve
                        && self.cliques[cid].num_zero_fixed
                            >= cpp_max(10, (self.cliques[cid].end - self.cliques[cid].start) >> 1)
                    {
                        let equality = self.cliques[cid].equality;
                        let origin = self.cliques[cid].origin;
                        let start = self.cliques[cid].start as usize;
                        let end = self.cliques[cid].end as usize;
                        clq.clear();
                        clq.extend_from_slice(&self.cliqueentries[start..end]);

                        // assert(computeNumDeleted(clq) == cliques[cliqueid].numZeroFixed)
                        // (lines 1668-1674) dropped under NDEBUG

                        self.remove_clique(*cliqueid, false);
                        let col_deleted = &self.col_deleted;
                        clq.retain(|x| {
                            !(globaldom
                                .is_fixed_to_val(x.col() as HighsInt, 1 - x.val() as HighsInt)
                                || col_deleted[x.col() as usize] != 0)
                        });
                        if clq.len() > 1 {
                            self.do_add_clique(&clq, equality, origin);
                        } else if equality {
                            if clq.is_empty()
                                || globaldom.is_fixed_to_val(
                                    clq[0].col() as HighsInt,
                                    1 - clq[0].val() as HighsInt,
                                )
                            {
                                globaldom.mark_infeasible(HighsDomain_Reason::unspecified());
                            } else if !globaldom.is_fixed(clq[0].col() as HighsInt) {
                                self.fix_col(ms, globaldom, clq[0].complement(), false);
                            }
                        }
                    }
                });
                continue;
            }

            if self.in_presolve {
                // during presolve we only count the number of zeros within each clique
                // and only remove fully redundant cliques that are larger than two
                // in the process since during presolve a lot of cliques of size two
                // may be found by probing and will be deleted upon rebuild anyways
                v_hash_lists.for_each(&mut |cliqueid: &HighsInt, _: &HighsInt| {
                    let cid = *cliqueid as usize;
                    self.cliques[cid].num_zero_fixed += 1;
                    if self.cliques[cid].num_active() <= 1 {
                        if x1b_fix_last {
                            self.fix_last_active_and_remove(ms, globaldom, *cliqueid);
                        } else {
                            self.remove_clique(*cliqueid, true);
                        }
                    }
                });
                continue;
            }

            v_hash_lists_size_two.for_each(&mut |cliqueid: &HighsInt, _: &()| {
                self.remove_clique(*cliqueid, true);
            });

            debug_assert!(self.cliquehitinds.is_empty());
            let mut clq: Vec<HighsCliqueTable_CliqueVar> = Vec::new();
            v_hash_lists.for_each(&mut |cliqueid: &HighsInt, _: &HighsInt| {
                let cid = *cliqueid as usize;
                // assert(cliqueentries[entry.value()].val == 1 - v.val);

                self.cliques[cid].num_zero_fixed += 1;
                if self.cliques[cid].num_active() <= 1 {
                    if x1b_fix_last {
                        self.fix_last_active_and_remove(ms, globaldom, *cliqueid);
                    } else {
                        self.remove_clique(*cliqueid, true);
                    }
                } else if self.cliques[cid].num_zero_fixed
                    >= cpp_max(10, (self.cliques[cid].end - self.cliques[cid].start) >> 1)
                {
                    let start = self.cliques[cid].start as usize;
                    let end = self.cliques[cid].end as usize;
                    clq.clear();
                    clq.extend_from_slice(&self.cliqueentries[start..end]);

                    // assert(computeNumDeleted(clq) == cliques[cliqueid].numZeroFixed);

                    self.remove_clique(*cliqueid, true);
                    clq.retain(|x| {
                        !(globaldom.is_fixed(x.col() as HighsInt)
                            && globaldom.col_lower_[x.col() as usize] == (1 - x.val()) as f64)
                    });
                    if clq.len() > 1 {
                        self.do_add_clique(&clq, false, kHighsIInf);
                    }
                }
            });
        }

        self.propagate_and_cleanup(ms, globaldom);
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::propagateAndCleanup (lines 1569-1593)
    /// C++: `void propagateAndCleanup(HighsDomain & globaldom)`
    /// PORT NOTE: signature changed to take `ms: &mut` (`propagate()` needs it, M1 ruling).
    pub fn propagate_and_cleanup(&mut self, ms: &mut HighsMipSolver, globaldom: &mut HighsDomain) {
        // The C++ holds a reference to the change stack; only its size is read, so the length is
        // re-read from the domain at each point.
        let mut start = globaldom.get_domain_change_stack().len();
        lend!(ms.mipdata_.cliquetable, self, globaldom.propagate(ms));
        let mut end = globaldom.get_domain_change_stack().len();

        while !globaldom.infeasible() && start != end {
            for k in start..end {
                let col = globaldom.get_domain_change_stack()[k].column;
                if !globaldom.is_fixed(col) {
                    continue;
                }
                if globaldom.col_lower_[col as usize] != 1.0
                    && globaldom.col_lower_[col as usize] != 0.0
                {
                    continue;
                }

                let fixval = globaldom.col_lower_[col as usize] as HighsInt;
                let v = HighsCliqueTable_CliqueVar::new(col, 1 - fixval);
                if self.num_cliques_2(v) != 0 {
                    self.vertex_infeasible(ms, globaldom, col, 1 - fixval);
                    if globaldom.infeasible() {
                        return;
                    }
                }
            }
            start = globaldom.get_domain_change_stack().len();
            lend!(ms.mipdata_.cliquetable, self, globaldom.propagate(ms));
            end = globaldom.get_domain_change_stack().len();
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::vertexInfeasible (lines 1595-1603)
    /// C++: `void vertexInfeasible(HighsDomain & globaldom, HighsInt col, HighsInt val)`
    /// PORT NOTE: signature changed to take `ms: &mut` (`fixCol` needs it, M1 ruling).
    pub fn vertex_infeasible(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldom: &mut HighsDomain,
        col: HighsInt,
        val: HighsInt,
    ) {
        let wasfixed = globaldom.is_fixed(col);
        lend!(
            ms.mipdata_.cliquetable,
            self,
            globaldom.fix_col(ms, col, (1 - val) as f64, HighsDomain_Reason::unspecified())
        );
        if globaldom.infeasible() {
            return;
        }
        if !wasfixed {
            self.nfixings += 1;
        }
        self.infeasvertexstack
            .push(HighsCliqueTable_CliqueVar::new(col, val));
        self.process_infeasible_vertices(ms, globaldom);
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::separateCliques (lines 1605-1714)
    /// C++: `void separateCliques(const HighsMipSolver & mipsolver, const std::vector<double> & sol, HighsCutPool & cutpool, double feastol, HighsRandom & randgen, int64_t & localNumNeighbourhoodQueries)`
    ///
    /// PORT NOTE: signature changed because (a) `cutpool.addCut` takes `ms: &mut` (cut_pool.rs),
    /// and (b) in the serial path HighsSeparation.cpp:110-116 passes the table's own `randgen`
    /// and `numNeighbourhoodQueries`, which Rust cannot alias with `&mut self`. `None` stands for
    /// the table's own member (`&randgen == &this->randgen` is `randgen.is_none()`); `Some(x)` is
    /// the worker's copy (parallel lock active).
    ///
    /// HiGHS quirk: `ADD_ZERO_WEIGHT_VARS` is `#define`d at HighsCliqueTable.cpp:25, so the
    /// zero-weight extension (`data.Z`) branch is the compiled one and is translated here.
    #[allow(clippy::too_many_arguments)]
    pub fn separate_cliques(
        &mut self,
        ms: &mut HighsMipSolver,
        sol: &[f64],
        cutpool: &mut HighsCutPool,
        feastol: f64,
        mut randgen: Option<&mut HighsRandom>,
        local_num_neighbourhood_queries: Option<&mut i64>,
    ) {
        let mut data = HighsCliqueTable_BronKerboschData::new(sol);
        data.feastol = feastol;
        data.max_neighbourhood_queries =
            1000000 + 100 * (ms.num_nonzero() as i64) + ms.mipdata_.total_lp_iterations * 1000;
        if self.num_neighbourhood_queries > data.max_neighbourhood_queries {
            return;
        }
        data.max_neighbourhood_queries -= self.num_neighbourhood_queries;

        for &i in ms.mipdata_.integral_cols.iter() {
            if self.colsubstituted[i as usize] != 0 || self.col_deleted[i as usize] != 0 {
                continue;
            }
            if self.num_cliques_3(i, false) != 0 {
                if HighsCliqueTable_CliqueVar::new(i, 0).weight(sol) > feastol {
                    data.p.push(HighsCliqueTable_CliqueVar::new(i, 0));
                } else {
                    data.z.push(HighsCliqueTable_CliqueVar::new(i, 0));
                }
            }
            if self.num_cliques_3(i, true) != 0 {
                if HighsCliqueTable_CliqueVar::new(i, 1).weight(sol) > feastol {
                    data.p.push(HighsCliqueTable_CliqueVar::new(i, 1));
                } else {
                    data.z.push(HighsCliqueTable_CliqueVar::new(i, 1));
                }
            }
        }

        let plen = data.p.len() as HighsInt;
        self.bron_kerbosch_recurse(&mut data, plen, &[]);

        let mut runcliquesubsumption = false;
        let mut inds: Vec<HighsInt> = Vec::new();
        let mut vals: Vec<f64> = Vec::new();
        for ci in 0..data.cliques.len() {
            let mut extensionend = data.z.len() as HighsInt;
            for vi in 0..data.cliques[ci].len() {
                let v = data.cliques[ci][vi];
                extensionend = self.partition_neighbourhood(
                    &mut data.neighbourhood_inds,
                    &mut data.num_neighbourhood_queries,
                    v,
                    &mut data.z[..extensionend as usize],
                );
                if extensionend == 0 {
                    break;
                }
            }

            if extensionend != 0 {
                match randgen.as_deref_mut() {
                    Some(r) => r.shuffle(&mut data.z[..], extensionend),
                    None => self.randgen.shuffle(&mut data.z[..], extensionend),
                }

                let mut i: HighsInt = 0;
                while i < extensionend {
                    let k = i + 1;
                    let zi = data.z[i as usize];
                    extensionend = k + self.partition_neighbourhood(
                        &mut data.neighbourhood_inds,
                        &mut data.num_neighbourhood_queries,
                        zi,
                        &mut data.z[k as usize..extensionend as usize],
                    );
                    i += 1;
                }

                let ext: Vec<HighsCliqueTable_CliqueVar> = data.z[..extensionend as usize].to_vec();
                data.cliques[ci].extend_from_slice(&ext);
            }

            let clique = &data.cliques[ci];
            let mut rhs: f64 = 1.0;
            runcliquesubsumption = self.cliques.len() > 2;
            inds.clear();
            vals.clear();

            for v in clique.iter() {
                inds.push(v.col() as HighsInt);
                if v.val() == 0 {
                    vals.push(-1.0);
                    rhs -= 1.0;
                } else {
                    vals.push(1.0);
                }
            }

            rhs = (rhs + 0.5).floor();

            let n = inds.len() as HighsInt;
            cutpool.add_cut(ms, &mut inds, &mut vals, n, rhs, true, false, false, false);
        }

        match local_num_neighbourhood_queries {
            Some(q) => *q += data.num_neighbourhood_queries,
            None => self.num_neighbourhood_queries += data.num_neighbourhood_queries,
        }

        if runcliquesubsumption && randgen.is_none() {
            for clique in data.cliques.iter_mut() {
                let nremoved = self.run_clique_subsumption(&ms.mipdata_.domains[0], clique);

                if clique.is_empty() {
                    continue;
                }
                if nremoved != 0 {
                    self.do_add_clique(&clique[..], false, -1);
                }
            }
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::separateCliques (lines 1716-1742)
    /// C++: `std::vector<std::vector<CliqueVar>> separateCliques(const std::vector<double> & sol, const HighsDomain & globaldom, double feastol)`
    pub fn separate_cliques_2(
        &mut self,
        sol: &[f64],
        globaldom: &HighsDomain,
        feastol: f64,
    ) -> Vec<Vec<HighsCliqueTable_CliqueVar>> {
        // PORT NOTE: the C++ body is `exit(0);` followed by an `#if 0` block (dropped). The
        // process exit is kept literally.
        std::process::exit(0)
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::computeMaximalCliques (lines 1744-1771)
    /// C++: `std::vector<std::vector<CliqueVar>> computeMaximalCliques(const std::vector<CliqueVar> & vars, double feastol)`
    pub fn compute_maximal_cliques(
        &mut self,
        vars: &[HighsCliqueTable_CliqueVar],
        feastol: f64,
    ) -> Vec<Vec<HighsCliqueTable_CliqueVar>> {
        // return if there are no variables
        if vars.is_empty() {
            return Vec::new();
        }

        // find max column index in clique variables
        let mut maxcolindex: usize = 0;
        for var in vars {
            maxcolindex = maxcolindex.max(var.col() as usize);
        }

        // set up data
        let mut sol: Vec<f64> = vec![0.0; maxcolindex + 1];
        for var in vars {
            sol[var.col() as usize] = var.val() as f64;
        }
        let mut data = HighsCliqueTable_BronKerboschData::new(&sol);
        data.feastol = feastol;

        for var in vars {
            if self.colsubstituted[var.col() as usize] != 0
                || self.col_deleted[var.col() as usize] != 0
            {
                continue;
            }
            if self.num_cliques_2(*var) != 0 {
                data.p.push(*var);
            }
        }

        let plen = data.p.len() as HighsInt;
        self.bron_kerbosch_recurse(&mut data, plen, &[]);

        data.cliques
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::addImplications (lines 1773-1827)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::addImplications (lines 1906-1962)
    /// C++: `void addImplications(HighsDomain & domain, HighsInt col, HighsInt val)`
    /// PORT NOTE: signature changed (re-entrancy):
    /// an associated fn reading `ms.mipdata_.cliquetable` afresh, because it is called from
    /// `HighsDomain::change_bound` and calls `change_bound` back (which re-enters this fn and
    /// `HighsImplications::apply_implications`). The table must be attached (or lent with
    /// `lend!`) when `change_bound` runs; nothing on this path mutates it.
    pub fn add_implications(
        ms: &mut HighsMipSolver,
        domain: &mut HighsDomain,
        col: HighsInt,
        val: HighsInt,
    ) {
        let mut v = HighsCliqueTable_CliqueVar::new(col, val);

        loop {
            let ct = &ms.mipdata_.cliquetable;
            if ct.colsubstituted[v.col() as usize] == 0 {
                break;
            }
            debug_assert!(
                ct.substitutions.len() as HighsInt > ct.colsubstituted[v.col() as usize] - 1
            );
            let subst =
                ct.substitutions[(ct.colsubstituted[v.col() as usize] - 1) as usize].clone();
            v = if v.val() == 1 {
                subst.replace
            } else {
                subst.replace.complement()
            };
            if v.val() == 1 {
                if domain.col_lower_[v.col() as usize] == 1.0 {
                    continue;
                }

                domain.change_bound(
                    ms,
                    HighsDomainChange {
                        boundval: 1.0,
                        column: v.col() as HighsInt,
                        boundtype: HighsBoundType::kLower,
                    },
                    HighsDomain_Reason::clique_table(col, val),
                );
                if domain.infeasible() {
                    return;
                }
            } else {
                if domain.col_upper_[v.col() as usize] == 0.0 {
                    continue;
                }

                domain.change_bound(
                    ms,
                    HighsDomainChange {
                        boundval: 0.0,
                        column: v.col() as HighsInt,
                        boundtype: HighsBoundType::kUpper,
                    },
                    HighsDomain_Reason::clique_table(col, val),
                );
                if domain.infeasible() {
                    return;
                }
            }
        }

        // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 1935-1936 at 7b890c255
        let clqp = ms.mipdata_.cliquetable.x1b_presolve_clique_table;
        // lambda `doFixings`, applied to every clique id of a list in turn until it
        // reports infeasibility (the `for_each` early exit).
        //
        // PORT NOTE (perf): the entries are scanned with the clique table and the
        // domain bounds borrowed once, across clique boundaries, instead of
        // re-reading them through `ms` for every clique and entry. Most cliques need
        // no fixing, so the scan usually runs through the whole list on one borrow.
        // Every `change_bound` ends the scan and the borrows are taken afresh after
        // it. Each clique's `start`/`end` are read once when its scan begins and the
        // entries are re-read after every `change_bound`, as in the C++ lambda, so
        // the entries tested, the calls made and their order are those of C++.
        //
        // The clique ids to visit are `add_implications_ids[k..kend]`.
        let do_fixings_list =
            |ms: &mut HighsMipSolver, domain: &mut HighsDomain, k: usize, kend: usize| -> bool {
                let mut k = k;
                let mut i = 0usize;
                let mut end = 0usize;
                loop {
                    // scan to the next entry that needs a bound change
                    let mut found: Option<HighsDomainChange> = None;
                    ADD_IMPLICATIONS_IDS.with(|cell| {
                        let ct = &ms.mipdata_.cliquetable;
                        let ids_guard = cell.borrow();
                        let ids = &ids_guard[..kend];
                        let cliques = &ct.cliques[..];
                        let entries = &ct.cliqueentries[..];
                        let col_deleted = &ct.col_deleted[..];
                        let col_lower = &domain.col_lower_[..];
                        let col_upper = &domain.col_upper_[..];
                        'scan: loop {
                            while i == end {
                                let Some(&cid) = ids.get(k) else {
                                    break 'scan;
                                };
                                k += 1;
                                i = cliques[cid as usize].start as usize;
                                end = cliques[cid as usize].end as usize;
                            }
                            for ce in &entries[i..end] {
                                i += 1;
                                if clqp {
                                    if ce.col() == v.col() || col_deleted[ce.col() as usize] != 0 {
                                        continue;
                                    }
                                } else if ce.col() == v.col() {
                                    continue;
                                }

                                if ce.val() == 1 {
                                    if col_upper[ce.col() as usize] == 0.0 {
                                        continue;
                                    }
                                    found = Some(HighsDomainChange {
                                        boundval: 0.0,
                                        column: ce.col() as HighsInt,
                                        boundtype: HighsBoundType::kUpper,
                                    });
                                } else {
                                    if col_lower[ce.col() as usize] == 1.0 {
                                        continue;
                                    }
                                    found = Some(HighsDomainChange {
                                        boundval: 1.0,
                                        column: ce.col() as HighsInt,
                                        boundtype: HighsBoundType::kLower,
                                    });
                                }
                                break 'scan;
                            }
                        }
                    });
                    let Some(boundchg) = found else {
                        return false;
                    };
                    domain.change_bound(ms, boundchg, HighsDomain_Reason::clique_table(col, val));
                    if domain.infeasible() {
                        return true;
                    }
                }
            };

        // PORT NOTE: C++ runs `doFixings` inside `invertedHashList[..].for_each`, i.e. while
        // the tree is borrowed. The ids are collected first (same traversal order) so the
        // scan can take `ms` mutably; the tree is not modified on this path, so the visit
        // sequence and the early exit are unchanged. The ids go on a stack shared by the
        // nested calls (`change_bound` re-enters this function when it fixes a binary):
        // each call pushes its ids above its caller's and truncates back to `base` when
        // done, so the buffer's allocation is reused at every depth.
        let base = ADD_IMPLICATIONS_IDS.with(|cell| cell.borrow().len());
        let kend = ADD_IMPLICATIONS_IDS.with(|cell| {
            let mut ids = cell.borrow_mut();
            ms.mipdata_.cliquetable.inverted_hash_list[v.index() as usize]
                .for_each(&mut |cid: &HighsInt, _: &HighsInt| ids.push(*cid));
            ids.len()
        });
        let infeas = do_fixings_list(ms, domain, base, kend);
        ADD_IMPLICATIONS_IDS.with(|cell| cell.borrow_mut().truncate(base));

        if !infeas {
            let kend = ADD_IMPLICATIONS_IDS.with(|cell| {
                let mut ids = cell.borrow_mut();
                ms.mipdata_.cliquetable.inverted_hash_list_size_two[v.index() as usize]
                    .for_each(&mut |cid: &HighsInt, _: &()| ids.push(*cid));
                ids.len()
            });
            do_fixings_list(ms, domain, base, kend);
            ADD_IMPLICATIONS_IDS.with(|cell| cell.borrow_mut().truncate(base));
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::cleanupFixed (lines 1829-1845)
    /// C++: `void cleanupFixed(HighsDomain & globaldom)`
    /// PORT NOTE: signature changed to take `ms: &mut` (M1 ruling).
    pub fn cleanup_fixed(&mut self, ms: &mut HighsMipSolver, globaldom: &mut HighsDomain) {
        let numcol = globaldom.col_upper_.len();
        let oldnfixings = self.nfixings;
        for i in 0..numcol {
            if self.col_deleted[i] != 0 || !globaldom.is_fixed(i as HighsInt) {
                continue;
            }
            if globaldom.col_lower_[i] != 1.0 && globaldom.col_lower_[i] != 0.0 {
                continue;
            }

            let fixval = globaldom.col_lower_[i] as HighsInt;
            let v = HighsCliqueTable_CliqueVar::new(i as HighsInt, 1 - fixval);

            self.vertex_infeasible(ms, globaldom, v.col() as HighsInt, v.val() as HighsInt);
            if globaldom.infeasible() {
                return;
            }
        }

        if self.nfixings != oldnfixings {
            self.propagate_and_cleanup(ms, globaldom);
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::getNumImplications (lines 1847-1864)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::getNumImplications (lines 1981-1998)
    /// C++: `HighsInt getNumImplications(HighsInt col) const`
    pub fn get_num_implications(&self, col: HighsInt) -> HighsInt {
        // first count all cliques as one implication, so that cliques of size two
        // are accounted for already
        let i0 = HighsCliqueTable_CliqueVar::new(col, 0).index();
        let i1 = HighsCliqueTable_CliqueVar::new(col, 1).index();
        let mut numimplics = self.numcliquesvar[i0 as usize] + self.numcliquesvar[i1 as usize];

        // now loop over cliques larger than size two and add the cliquelength - 1 as
        // additional implications
        // PORT NOTE: the Rust hash tree passes (key, value); the C++ one-argument lambda
        // reads only the key (`cliqueid`), which is what `cliqueid` is here.
        // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 1991, 2009 at 7b890c255
        let clqp = self.x1b_presolve_clique_table;
        let mut count_implics = |cliqueid: &HighsInt, _value: &HighsInt| {
            let mut nimplics = if clqp {
                self.cliques[*cliqueid as usize].num_active() - 1
            } else {
                self.cliques[*cliqueid as usize].end - self.cliques[*cliqueid as usize].start - 1
            };
            nimplics *= 1 + self.cliques[*cliqueid as usize].equality as HighsInt;
            numimplics += nimplics - 1;
        };
        self.inverted_hash_list[i0 as usize].for_each(&mut count_implics);
        self.inverted_hash_list[i1 as usize].for_each(&mut count_implics);
        numimplics
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::getNumImplications (lines 1866-1881)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::getNumImplications (lines 2000-2015)
    /// C++: `HighsInt getNumImplications(HighsInt col, bool val) const`
    pub fn get_num_implications_2(&self, col: HighsInt, val: bool) -> HighsInt {
        let i_val = HighsCliqueTable_CliqueVar::new(col, val as HighsInt).index();

        // each size two clique is one implication
        let mut numimplics = self.numcliquesvar[i_val as usize];

        // now loop over cliques larger than size two and add the cliquelength - 1 as
        // additional implications
        // PORT NOTE: see `get_num_implications` (the Rust hash tree passes (key, value)).
        // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 1991, 2009 at 7b890c255
        let clqp = self.x1b_presolve_clique_table;
        let mut count_implics = |cliqueid: &HighsInt, _value: &HighsInt| {
            let mut nimplics = if clqp {
                self.cliques[*cliqueid as usize].num_active() - 1
            } else {
                self.cliques[*cliqueid as usize].end - self.cliques[*cliqueid as usize].start - 1
            };
            nimplics *= 1 + self.cliques[*cliqueid as usize].equality as HighsInt;
            numimplics += nimplics - 1;
        };
        self.inverted_hash_list[i_val as usize].for_each(&mut count_implics);
        numimplics
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::runCliqueMerging (lines 1883-1985)
    /// C++: `void runCliqueMerging(HighsDomain & globaldomain, std::vector<CliqueVar> & clique, bool equation = false)`
    /// C++ default arguments: `equation = false`
    /// PORT NOTE: signature changed to take `ms: &mut` (M1 ruling). The C++ member
    /// `numNeighbourhoodQueries` is copied out around each `shrinkToNeighbourhood` call, since the
    /// callee cannot borrow `self` twice.
    pub fn run_clique_merging(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldomain: &mut HighsDomain,
        clique: &mut Vec<HighsCliqueTable_CliqueVar>,
        equation: bool,
    ) {
        let mut extensionstart = HighsCliqueTable_CliqueVar::default();
        let mut numcliques = kHighsIInf;
        self.iscandidate.resize(self.inverted_hash_list.len(), 0);
        let mut neighbourhood_inds: Vec<HighsInt> =
            Vec::with_capacity(self.inverted_hash_list.len());

        let initial_clique_size = clique.len();
        for i in 0..initial_clique_size {
            if globaldomain.is_fixed(clique[i].col() as HighsInt) {
                continue;
            }

            let this_num_clqs = self.num_cliques_2(clique[i]);
            if this_num_clqs < numcliques {
                numcliques = this_num_clqs;
                extensionstart = clique[i];
            }
        }

        if numcliques == kHighsIInf {
            clique.clear();
            return;
        }

        for i in 0..initial_clique_size {
            self.iscandidate[clique[i].index() as usize] = 1;
        }

        let mut add_cands = |cliqueid: HighsInt| {
            let start = self.cliques[cliqueid as usize].start;
            let end = self.cliques[cliqueid as usize].end;

            for i in start..end {
                let ce = self.cliqueentries[i as usize];
                if self.iscandidate[ce.index() as usize] != 0
                    || globaldomain.is_fixed(ce.col() as HighsInt)
                {
                    continue;
                }

                self.iscandidate[ce.index() as usize] = 1;
                clique.push(ce);
            }
        };

        self.inverted_hash_list[extensionstart.index() as usize]
            .for_each(&mut |cid: &HighsInt, _: &HighsInt| add_cands(*cid));

        self.inverted_hash_list_size_two[extensionstart.index() as usize]
            .for_each(&mut |cid: &HighsInt, _: &()| add_cands(*cid));

        let size_with_candidates = clique.len();
        for i in 0..size_with_candidates {
            self.iscandidate[clique[i].index() as usize] = 0;
        }

        let mut i = 0usize;
        while i != initial_clique_size && initial_clique_size < clique.len() {
            if clique[i] == extensionstart {
                i += 1;
                continue;
            }

            let v = clique[i];
            let mut nq = self.num_neighbourhood_queries;
            let new_size = initial_clique_size
                + self.shrink_to_neighbourhood(
                    &mut neighbourhood_inds,
                    &mut nq,
                    v,
                    &mut clique[initial_clique_size..],
                ) as usize;
            self.num_neighbourhood_queries = nq;
            clique.truncate(new_size);
            i += 1;
        }

        if initial_clique_size < clique.len() {
            // todo, shuffle extension vars?
            let n = (clique.len() - initial_clique_size) as HighsInt;
            self.randgen.shuffle(&mut clique[initial_clique_size..], n);
            let mut i = initial_clique_size;
            while i < clique.len() {
                let extvar = clique[i];
                i += 1;

                let mut nq = self.num_neighbourhood_queries;
                let new_size = i + self.shrink_to_neighbourhood(
                    &mut neighbourhood_inds,
                    &mut nq,
                    extvar,
                    &mut clique[i..],
                ) as usize;
                self.num_neighbourhood_queries = nq;
                clique.truncate(new_size);
            }
        }

        if equation {
            for i in initial_clique_size..clique.len() {
                let cv = clique[i];
                self.vertex_infeasible(
                    ms,
                    globaldomain,
                    cv.col() as HighsInt,
                    cv.val() as HighsInt,
                );
            }
        } else {
            self.run_clique_subsumption(globaldomain, clique);

            if !clique.is_empty() {
                clique.retain(|v| {
                    !(globaldomain.is_fixed(v.col() as HighsInt)
                        && (globaldomain.col_lower_[v.col() as usize] as i32)
                            == ((1 - v.val()) as i32))
                });
            }
        }

        self.process_infeasible_vertices(ms, globaldomain);
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::runCliqueMerging (lines 1987-2161)
    /// C++: `void runCliqueMerging(HighsDomain & globaldomain)`
    ///
    /// PORT NOTE: signature changed because `vertexInfeasible` and `processInfeasibleVertices`
    /// take `ms` (M1).
    ///
    /// PORT NOTE: the `removeCliques` callback of `cliqueSubsumption` would remove cliques while
    /// `cliqueSubsumption` holds `&mut self`. The callback records the ids and they are removed
    /// right after the call, in the same order. This is equivalent: every clique id occurs once in
    /// `cliquehitinds`, so the loop never reads a clique after its callback; `removeClique(x)`
    /// only erases `x` from the hash lists (the loop's `contains(cliqueid)` asks about other ids)
    /// and appends to `deletedrows`/`freeslots`/`freespaces`, which the loop does not read; the
    /// `infeasvertexstack` pushes of the loop go to a different vector.
    pub fn run_clique_merging_2(
        &mut self,
        ms: &mut HighsMipSolver,
        globaldomain: &mut HighsDomain,
    ) {
        let mut extensionvars: Vec<HighsCliqueTable_CliqueVar> = Vec::new();
        self.iscandidate.resize(self.inverted_hash_list.len(), 0);
        let mut neighbourhood_inds: Vec<HighsInt> =
            Vec::with_capacity(self.inverted_hash_list.len());

        let numcliqueslots = self.cliques.len() as HighsInt;
        let max_new_entries = self.num_entries + globaldomain.num_model_nonzeros(ms);
        let mut have_non_model_cliques_to_merge = false;
        for k in 0..numcliqueslots {
            let ku = k as usize;
            if self.cliques[ku].start == -1 {
                continue;
            }
            if !self.cliques[ku].equality && self.cliques[ku].origin == kHighsIInf {
                continue;
            }
            if self.cliques[ku].origin == -1 {
                have_non_model_cliques_to_merge = true;
                continue;
            }
            let numclqvars = self.cliques[ku].end - self.cliques[ku].start;
            debug_assert!(numclqvars != 0);
            if numclqvars == 0 {
                continue;
            }

            // C++ `CliqueVar* clqvars = &cliqueentries[cliques[k].start]`; `cliqueentries` is
            // not written until the clique is removed below, so `clqvars[i]` is read by index.
            let clqstart = self.cliques[ku].start as usize;
            let numclqvars = numclqvars as usize;

            let mut extensionstart = self.cliqueentries[clqstart];
            let mut numcliques = self.num_cliques_2(self.cliqueentries[clqstart]);
            for i in 1..numclqvars {
                let this_num_clqs = self.num_cliques_2(self.cliqueentries[clqstart + i]);
                if this_num_clqs < numcliques {
                    numcliques = this_num_clqs;
                    extensionstart = self.cliqueentries[clqstart + i];
                }
            }

            for i in 0..numclqvars {
                let idx = self.cliqueentries[clqstart + i].index() as usize;
                self.iscandidate[idx] = 1;
            }

            let mut add_cands = |cliqueid: HighsInt| {
                let start = self.cliques[cliqueid as usize].start;
                let end = self.cliques[cliqueid as usize].end;

                for i in start..end {
                    let ce = self.cliqueentries[i as usize];
                    if self.iscandidate[ce.index() as usize] != 0
                        || globaldomain.is_fixed(ce.col() as HighsInt)
                    {
                        continue;
                    }

                    self.iscandidate[ce.index() as usize] = 1;
                    extensionvars.push(ce);
                }
            };

            self.inverted_hash_list[extensionstart.index() as usize]
                .for_each(&mut |cid: &HighsInt, _: &HighsInt| add_cands(*cid));

            self.inverted_hash_list_size_two[extensionstart.index() as usize]
                .for_each(&mut |cid: &HighsInt, _: &()| add_cands(*cid));

            for i in 0..numclqvars {
                let idx = self.cliqueentries[clqstart + i].index() as usize;
                self.iscandidate[idx] = 0;
            }
            for v in extensionvars.iter() {
                self.iscandidate[v.index() as usize] = 0;
            }

            let mut i = 0usize;
            while i != numclqvars && !extensionvars.is_empty() {
                let clqvar = self.cliqueentries[clqstart + i];
                if clqvar == extensionstart {
                    i += 1;
                    continue;
                }

                let mut nq = self.num_neighbourhood_queries;
                let new_size = self.shrink_to_neighbourhood(
                    &mut neighbourhood_inds,
                    &mut nq,
                    clqvar,
                    &mut extensionvars[..],
                );
                self.num_neighbourhood_queries = nq;
                extensionvars.truncate(new_size as usize);
                i += 1;
            }

            if !extensionvars.is_empty() {
                // todo, shuffle extension vars?
                let n = extensionvars.len() as HighsInt;
                self.randgen.shuffle(&mut extensionvars[..], n);
                let mut i = 0usize;
                while i < extensionvars.len() {
                    let extvar = extensionvars[i];
                    i += 1;

                    let mut nq = self.num_neighbourhood_queries;
                    let new_size = i + self.shrink_to_neighbourhood(
                        &mut neighbourhood_inds,
                        &mut nq,
                        extvar,
                        &mut extensionvars[i..],
                    ) as usize;
                    self.num_neighbourhood_queries = nq;
                    extensionvars.truncate(new_size);
                }
            }

            if self.cliques[ku].equality {
                for idx in 0..extensionvars.len() {
                    let v = extensionvars[idx];
                    self.vertex_infeasible(
                        ms,
                        globaldomain,
                        v.col() as HighsInt,
                        v.val() as HighsInt,
                    );
                }
            } else {
                let originrow = self.cliques[ku].origin;
                self.cliques[ku].origin = kHighsIInf;

                let num_extensions = extensionvars.len();
                let cstart = self.cliques[ku].start as usize;
                let cend = self.cliques[ku].end as usize;
                extensionvars.extend_from_slice(&self.cliqueentries[cstart..cend]);
                // C++ erase(remove_if(begin() + numExtensions, end(), colDeleted)): a stable
                // compaction of the appended tail.
                let tail = extensionvars.split_off(num_extensions);
                extensionvars.extend(
                    tail.into_iter()
                        .filter(|clqvar| self.col_deleted[clqvar.col() as usize] == 0),
                );
                self.remove_clique(k, true);

                // callback for removing cliques (deferred, see PORT NOTE)
                let mut to_remove: Vec<HighsInt> = Vec::new();

                // check extension
                let mut redundant = false;
                let mut dominating_origin: HighsInt = kHighsIInf;
                self.clique_subsumption(
                    &extensionvars,
                    &mut redundant,
                    &mut dominating_origin,
                    &mut |cliqueid| to_remove.push(cliqueid),
                );
                for cliqueid in to_remove {
                    self.remove_clique(cliqueid, true);
                }

                if !redundant {
                    for i in 0..num_extensions {
                        self.cliqueextensions.push((originrow, extensionvars[i]));
                    }

                    extensionvars.retain(|v| {
                        !(globaldomain.is_fixed(v.col() as HighsInt)
                            && (globaldomain.col_lower_[v.col() as usize] as i32)
                                == (1 - v.val() as i32))
                    });

                    if extensionvars.len() > 1 {
                        self.do_add_clique(&extensionvars, false, originrow);
                    }
                } else {
                    // the extended clique is redundant, check if the row can be removed
                    if dominating_origin != kHighsIInf {
                        self.deletedrows.push(originrow);
                    } else {
                        // this clique is redundant in the cliquetable but its row is not
                        // necessarily. Also there might be rows that have been deleted due to
                        // being dominated by this row after adding the lifted entries so they
                        // must be added to the cliqueextension vector
                        for i in 0..num_extensions {
                            self.cliqueextensions.push((originrow, extensionvars[i]));
                        }
                    }
                }
            }

            extensionvars.clear();
            self.process_infeasible_vertices(ms, globaldomain);

            if self.num_entries >= max_new_entries {
                break;
            }
            // printf("nonzeroDelta: %d, maxNonzeroDelta: %d\n", nonzeroDelta,
            // maxNonzeroDelta);
        }

        if have_non_model_cliques_to_merge {
            for k in 0..numcliqueslots {
                let ku = k as usize;
                if self.cliques[ku].start == -1 {
                    continue;
                }
                if self.cliques[ku].origin != -1 {
                    continue;
                }
                // if (cliques[k].end - cliques[k].start <= 1000) continue;

                // printf("numEntries before: %d\n", numEntries);
                extensionvars.clear();
                let cstart = self.cliques[ku].start as usize;
                let cend = self.cliques[ku].end as usize;
                extensionvars.extend_from_slice(&self.cliqueentries[cstart..cend]);
                self.remove_clique(k, true);
                self.run_clique_merging(ms, globaldomain, &mut extensionvars, false);

                if extensionvars.len() > 1 {
                    self.do_add_clique(&extensionvars, false, kHighsIInf);
                }
            }
        }
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::rebuild (lines 2163-2204)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::rebuild (lines 2291-2331)
    /// C++: `void rebuild(HighsInt ncols, const presolve::HighsPostsolveStack & postSolveStack, const HighsDomain & globaldomain, const std::vector<HighsInt> & orig2reducedcol, const std::vector<HighsInt> & orig2reducedrow)`
    /// PORT NOTE: signature gained `ms: &mut HighsMipSolver` (for `isBinary`, M1 ruling; `isBinary`
    /// takes `&mut` in domain/mod.rs:592).
    pub fn rebuild(
        &mut self,
        ms: &mut HighsMipSolver,
        ncols: HighsInt,
        post_solve_stack: &HighsPostsolveStack,
        globaldomain: &HighsDomain,
        orig2reducedcol: &[HighsInt],
        orig2reducedrow: &[HighsInt],
    ) {
        let mut new_clique_table = HighsCliqueTable::new(ncols);
        new_clique_table.set_presolve_flag(self.in_presolve);
        // HiGHS 1.16 #3540 (h16_clique_links): plumbing, flag copied to the rebuilt table
        new_clique_table.copy_h16_flags_from(self);
        new_clique_table.set_min_entries_for_parallelism(self.min_entries_for_parallelism);
        // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 2300-2322 at 7b890c255
        let clqp = self.x1b_presolve_clique_table;
        for i in 0..self.cliques.len() {
            if self.cliques[i].start == -1 {
                continue;
            }
            let start = self.cliques[i].start as usize;
            let end = self.cliques[i].end as usize;
            let oldnumvars = (end - start) as HighsInt;
            // C++ latest: `const HighsInt numActive = cliques[i].numActive();`
            let num_active = if clqp {
                self.cliques[i].num_active()
            } else {
                0
            };

            for k in start..end {
                let col = orig2reducedcol[self.cliqueentries[k].col() as usize];

                if clqp {
                    if self.col_deleted[self.cliqueentries[k].col() as usize] != 0
                        || col == -1
                        || !globaldomain.is_binary(ms, col)
                        || !post_solve_stack.is_col_linearly_transformable(col)
                    {
                        self.cliqueentries[k].set_col(kHighsIInf as HighsUInt);
                    } else {
                        self.cliqueentries[k].set_col(col as HighsUInt);
                    }
                } else if col == -1
                    || !globaldomain.is_binary(ms, col)
                    || !post_solve_stack.is_col_linearly_transformable(col)
                {
                    self.cliqueentries[k].set_col(kHighsIInf as HighsUInt);
                } else {
                    self.cliqueentries[k].set_col(col as HighsUInt);
                }
            }

            // std::remove_if(cliqueentries[start..end], col == kHighsIInf), order preserving
            let mut w = start;
            for k in start..end {
                if self.cliqueentries[k].col() != kHighsIInf as HighsUInt {
                    self.cliqueentries[w] = self.cliqueentries[k];
                    w += 1;
                }
            }
            let numvars = (w - start) as HighsInt;

            // since we do not know how variables in the clique that have been deleted
            // are replaced (i.e. are they fixed to 0 or 1, or substituted) we relax
            // them out which means the equality status needs to be set to false
            if numvars <= 1 {
                continue;
            }

            let origin = if self.cliques[i].origin != kHighsIInf {
                -1
            } else {
                kHighsIInf
            };
            let equality = if clqp {
                // Preserve equality if all active entries survive, else relax to <=.
                numvars == num_active && self.cliques[i].equality
            } else if numvars != oldnumvars {
                false
            } else {
                self.cliques[i].equality
            };
            new_clique_table.do_add_clique(
                &self.cliqueentries[start..start + numvars as usize],
                equality,
                origin,
            );
        }

        new_clique_table.set_allow_parallel(self.allow_parallel);
        *self = new_clique_table;
    }

    /// HiGHS: HighsCliqueTable.cpp::HighsCliqueTable::buildFrom (lines 2206-2242)
    /// HiGHS latest (#3312, 7b890c255): HighsCliqueTable.cpp::HighsCliqueTable::buildFrom (lines 2333-2370)
    /// C++: `void buildFrom(const HighsLp * origModel, const HighsCliqueTable & init)`
    pub fn build_from(&mut self, orig_model: Option<&HighsLp>, init: &HighsCliqueTable) {
        debug_assert!(init.colsubstituted.len() == self.colsubstituted.len());
        let ncols = init.colsubstituted.len() as HighsInt;
        let mut new_clique_table = HighsCliqueTable::new(ncols);
        new_clique_table.set_presolve_flag(self.in_presolve);
        // HiGHS 1.16 #3540 (h16_clique_links): plumbing, flag copied to the new table
        new_clique_table.copy_h16_flags_from(self);
        new_clique_table.set_min_entries_for_parallelism(self.min_entries_for_parallelism);
        // PORT NOTE: C++ dereferences `origModel` unconditionally.
        let orig_model = orig_model.unwrap();
        // HiGHS #3312 (mip_presolve_clique_table): the gate, read from this table (as the
        // HiGHS 1.16 flags are copied from it)
        let clqp = self.x1b_presolve_clique_table;
        let ncliques = init.cliques.len();
        let mut clq_buffer: Vec<HighsCliqueTable_CliqueVar> =
            Vec::with_capacity((2 * orig_model.num_col_) as usize);
        for i in 0..ncliques {
            if init.cliques[i].start == -1 {
                continue;
            }
            if init.cliques[i].num_active() <= 1 {
                continue;
            }

            let start = init.cliques[i].start as usize;
            let end = init.cliques[i].end as usize;
            clq_buffer.clear();
            clq_buffer.extend_from_slice(&init.cliqueentries[start..end]);
            if clqp {
                // HiGHS #3312 (mip_presolve_clique_table): HighsCliqueTable.cpp lines 2349-2356 at 7b890c255
                clq_buffer.retain(|v| {
                    !(orig_model.col_lower_[v.col() as usize] != 0.0
                        || orig_model.col_upper_[v.col() as usize] != 1.0
                        || init.col_deleted[v.col() as usize] != 0)
                });
            } else {
                clq_buffer.retain(|v| {
                    !(orig_model.col_lower_[v.col() as usize] != 0.0
                        || orig_model.col_upper_[v.col() as usize] != 1.0)
                });
            }
            if clq_buffer.len() <= 1 {
                continue;
            }

            let origin = if init.cliques[i].origin != kHighsIInf {
                -1
            } else {
                kHighsIInf
            };
            new_clique_table.do_add_clique(&clq_buffer, false, origin);
        }

        new_clique_table.colsubstituted = init.colsubstituted.clone();
        new_clique_table.substitutions = init.substitutions.clone();
        // Currently assume buildFrom is always used for sub-mips
        new_clique_table.set_allow_parallel(false);
        *self = new_clique_table;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_ms_clique5() -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 5;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0; 5];
        lp.col_lower_ = vec![0.0; 5];
        lp.col_upper_ = vec![1.0; 5];
        lp.integrality_ = vec![HighsVarType::kInteger; 5];
        lp.a_matrix_.num_col_ = 5;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0; 6];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.mipdata_.post_solve_stack.initialize_index_maps(0, 5);
        ms.mipdata_.cliquetable = Slot::new(HighsCliqueTable::new(5));
        let imp = crate::mip::implications::HighsImplications::new(&ms);
        *ms.mipdata_.implications = imp;
        ms.mipdata_.domains = vec![Slot::new(HighsDomain::new(&ms))];
        ms
    }

    /// Port of the upstream test "MIP-equality-clique-fixing" (and its val=0
    /// variant): an equality clique over five columns, four of them fixed
    /// through vertexInfeasible. Returns the domain after the last fixing.
    fn run_x1b_flar_clique(
        fix_last: bool,
        in_presolve: bool,
        clique_val: HighsInt,
        fixed_val: f64,
    ) -> HighsDomain {
        let mut ms = build_ms_clique5();
        ms.options_mip_.mip_fix_last_active = fix_last;
        let mut ct = HighsCliqueTable::new(5);
        let clique: Vec<HighsCliqueTable_CliqueVar> = (0..5)
            .map(|i| HighsCliqueTable_CliqueVar::new(i, clique_val))
            .collect();
        ct.do_add_clique(&clique, true, kHighsIInf);
        // upstream test: cliquetable.setPresolveFlag(true)
        ct.set_presolve_flag(in_presolve);
        let mut dom = HighsDomain::new(&ms);
        for i in 0..4 {
            dom.fix_col(&mut ms, i, fixed_val, HighsDomain_Reason::unspecified());
            ct.vertex_infeasible(&mut ms, &mut dom, i, clique_val);
        }
        dom
    }

    #[test]
    fn x1b_flar_equality_clique_fixing() {
        // in_presolve = true is the upstream test; false covers the second call site
        for in_presolve in [true, false] {
            // on: the last active variable x4 is fixed to 1
            let dom_on = run_x1b_flar_clique(true, in_presolve, 1, 0.0);
            assert!(dom_on.is_fixed(4));
            assert_eq!(dom_on.col_lower_[4], 1.0);

            // off: 1.15.1 result, x4 is not fixed
            let dom_off = run_x1b_flar_clique(false, in_presolve, 1, 0.0);
            assert!(!dom_off.is_fixed(4));
        }
    }

    #[test]
    fn x1b_flar_equality_clique_fixing_to_zero() {
        // in_presolve = true is the upstream test; false covers the second call site
        for in_presolve in [true, false] {
            // on: the last active variable x4 is fixed to 0
            let dom_on = run_x1b_flar_clique(true, in_presolve, 0, 1.0);
            assert!(dom_on.is_fixed(4));
            assert_eq!(dom_on.col_upper_[4], 0.0);

            // off: 1.15.1 result, x4 is not fixed
            let dom_off = run_x1b_flar_clique(false, in_presolve, 0, 1.0);
            assert!(!dom_off.is_fixed(4));
        }
    }

    /// Port of 7b890c255:check/TestMipSolver.cpp:1614-1706 (no setPresolveFlag call: #3312
    /// removed it from the tests). `mip_fix_last_active` stays off, so the on result comes
    /// from the #3312 branch of `process_infeasible_vertices` alone.
    fn run_x1b_clqp_clique(
        clqp: bool,
        in_presolve_probing: bool,
        clique_val: HighsInt,
        fixed_val: f64,
    ) -> HighsDomain {
        let mut ms = build_ms_clique5();
        ms.options_mip_.mip_fix_last_active = false;
        ms.options_mip_.mip_presolve_clique_table = clqp;
        let mut ct = HighsCliqueTable::new(5);
        ct.set_h16_flags(&ms.options_mip_);
        let clique: Vec<HighsCliqueTable_CliqueVar> = (0..5)
            .map(|i| HighsCliqueTable_CliqueVar::new(i, clique_val))
            .collect();
        ct.do_add_clique(&clique, true, kHighsIInf);
        ct.set_presolve_flag(in_presolve_probing);
        let mut dom = HighsDomain::new(&ms);
        for i in 0..4 {
            dom.fix_col(&mut ms, i, fixed_val, HighsDomain_Reason::unspecified());
            ct.vertex_infeasible(&mut ms, &mut dom, i, clique_val);
        }
        dom
    }

    #[test]
    fn x1b_clqp_b_equality_clique_fixing() {
        // false is the upstream test (inPresolveProbing default); true covers probing
        for in_presolve_probing in [false, true] {
            // on: the last active variable x4 is fixed to 1
            let dom_on = run_x1b_clqp_clique(true, in_presolve_probing, 1, 0.0);
            assert!(!dom_on.infeasible());
            assert!(dom_on.is_fixed(4));
            assert_eq!(dom_on.col_lower_[4], 1.0);

            // off: 1.15.1 result, x4 is not fixed
            let dom_off = run_x1b_clqp_clique(false, in_presolve_probing, 1, 0.0);
            assert!(!dom_off.is_fixed(4));
        }
    }

    #[test]
    fn x1b_clqp_b_equality_clique_fixing_to_zero() {
        for in_presolve_probing in [false, true] {
            // on: the last active variable x4 is fixed to 0
            let dom_on = run_x1b_clqp_clique(true, in_presolve_probing, 0, 1.0);
            assert!(!dom_on.infeasible());
            assert!(dom_on.is_fixed(4));
            assert_eq!(dom_on.col_upper_[4], 0.0);

            // off: 1.15.1 result, x4 is not fixed
            let dom_off = run_x1b_clqp_clique(false, in_presolve_probing, 0, 1.0);
            assert!(!dom_off.is_fixed(4));
        }
    }

    #[test]
    fn x1b_clqp_b_remove_clique_record_deleted_row() {
        // removeClique(id, false) (7b890c255 HighsCliqueTable.cpp:818-821) does not push the
        // clique's origin row; the default `true` does, as in 1.15.1.
        for record in [true, false] {
            let mut ct = HighsCliqueTable::new(3);
            let clique: Vec<HighsCliqueTable_CliqueVar> = (0..3)
                .map(|i| HighsCliqueTable_CliqueVar::new(i, 1))
                .collect();
            ct.do_add_clique(&clique, false, 7);
            ct.remove_clique(0, record);
            assert_eq!(ct.cliques[0].start, -1);
            if record {
                assert_eq!(ct.deletedrows, vec![7]);
            } else {
                assert!(ct.deletedrows.is_empty());
            }
        }
    }

    fn build_ms(row_upper: f64, row_lower: f64, vals: Vec<f64>) -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0];
        lp.row_lower_ = vec![row_lower];
        lp.row_upper_ = vec![row_upper];
        lp.integrality_ = vec![HighsVarType::kInteger; 2];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 1];
        lp.a_matrix_.value_ = vals.clone();
        ms.mipdata_.a_rstart_ = vec![0, 2];
        ms.mipdata_.a_rindex_ = vec![0, 1];
        ms.mipdata_.a_rvalue_ = vals;
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.mipdata_.post_solve_stack.initialize_index_maps(1, 2);
        ms.mipdata_.domains = vec![Slot::new(HighsDomain::new(&ms))];
        ms
    }

    #[test]
    fn x1_clqn_a_extract_complemented_setppc() {
        // x0 - x1 <= 0: rhs 0 == 1 - numComp, clique {x0 = 1, x1 = 0}. With
        // transform_rows = false only the setppc block can extract it.
        let mut ms = build_ms(0.0, -kHighsInf, vec![1.0, -1.0]);
        ms.options_mip_.mip_clique_normalise = true;
        let mut ct = HighsCliqueTable::new(2);
        ct.extract_cliques_2(&mut ms, false);
        assert_eq!(ct.num_cliques_2(HighsCliqueTable_CliqueVar::new(0, 1)), 1);
        assert_eq!(ct.num_cliques_2(HighsCliqueTable_CliqueVar::new(1, 0)), 1);

        let mut ms_off = build_ms(0.0, -kHighsInf, vec![1.0, -1.0]);
        ms_off.options_mip_.mip_clique_normalise = false;
        let mut ct_off = HighsCliqueTable::new(2);
        ct_off.extract_cliques_2(&mut ms_off, false);
        assert_eq!(
            ct_off.num_cliques_2(HighsCliqueTable_CliqueVar::new(0, 1)),
            0
        );
        assert_eq!(
            ct_off.num_cliques_2(HighsCliqueTable_CliqueVar::new(1, 0)),
            0
        );
    }

    #[test]
    fn x1_clqn_a_plain_setppc_unchanged() {
        let mut ms_on = build_ms(1.0, -kHighsInf, vec![1.0, 1.0]);
        ms_on.options_mip_.mip_clique_normalise = true;
        let mut ct_on = HighsCliqueTable::new(2);
        ct_on.extract_cliques_2(&mut ms_on, true);

        let mut ms_off = build_ms(1.0, -kHighsInf, vec![1.0, 1.0]);
        ms_off.options_mip_.mip_clique_normalise = false;
        let mut ct_off = HighsCliqueTable::new(2);
        ct_off.extract_cliques_2(&mut ms_off, true);

        assert_eq!(
            ct_on.num_cliques_2(HighsCliqueTable_CliqueVar::new(0, 1)),
            ct_off.num_cliques_2(HighsCliqueTable_CliqueVar::new(0, 1))
        );
        assert_eq!(
            ct_on.num_cliques_2(HighsCliqueTable_CliqueVar::new(1, 1)),
            ct_off.num_cliques_2(HighsCliqueTable_CliqueVar::new(1, 1))
        );
    }
}
