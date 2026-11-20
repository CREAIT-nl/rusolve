//! HiGHS source: extern/amd/amd_2.c (lines 1281-1526)
//!
//! Translated from SuiteSparse AMD as vendored by HiGHS in `extern/amd` (BSD-3-Clause,
//! Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff).
//! See THIRD_PARTY_NOTICES.md.
//!
//! SUPERVARIABLE DETECTION (hash buckets on `hash % n`), RESTORE DEGREE LISTS
//! AND REMOVE NONPRINCIPAL SUPERVARIABLES FROM ELEMENT, and FINALIZE THE NEW
//! ELEMENT (statistics in `double`, kept in the C order).

use super::super::EMPTY;
use super::super::FLIP;
use super::super::{amd_int, amd_uint};
use super::Amd2;

impl<'a> Amd2<'a> {
    /// HiGHS: amd_2.c::Highs_amd_2 SUPERVARIABLE DETECTION + RESTORE + FINALIZE (lines 1281-1526)
    pub fn superv(&mut self) {
        // SUPERVARIABLE DETECTION (lines 1281-1419)
        self.pme = self.pme1;
        while self.pme <= self.pme2 {
            self.i = self.Iw[self.pme as usize];

            if self.Nv[self.i as usize] < 0 {
                // i is a principal variable in Lme

                // examine all hash buckets with 2 or more variables.  We do
                // this by examing all unique hash keys for supervariables in
                // the pattern Lme of the current element, me

                // let i = head of hash bucket, and empty the hash bucket

                self.hash = self.Last[self.i as usize] as amd_uint;

                // if Hhead array is not used:
                self.j = self.Head[self.hash as usize];
                if self.j == EMPTY {
                    // hash bucket and degree list are both empty
                    self.i = EMPTY;
                } else if self.j < EMPTY {
                    // degree list is empty
                    self.i = FLIP(self.j);
                    self.Head[self.hash as usize] = EMPTY;
                } else {
                    // degree list is not empty, restore Last [j] of head j
                    self.i = self.Last[self.j as usize];
                    self.Last[self.j as usize] = EMPTY;
                }

                // if separate Hhead array is used: (lines 1324-1327, commented out)

                while self.i != EMPTY && self.Next[self.i as usize] != EMPTY {
                    // this bucket has one or more variables following i.
                    // scan all of them to see if i can absorb any entries
                    // that follow i in hash bucket.  Scatter i into w.

                    self.ln = self.Len[self.i as usize];
                    self.eln = self.Elen[self.i as usize];

                    // do not flag the first element in the list (me)
                    self.p = self.Pe[self.i as usize] + 1;
                    while self.p < self.Pe[self.i as usize] + self.ln {
                        // PORT NOTE: `p <= x + ln - 1`, integers, clippy int_plus_one
                        self.W[self.Iw[self.p as usize] as usize] = self.wflg;
                        self.p += 1;
                    }

                    // scan every other entry j following i in bucket

                    self.jlast = self.i;
                    self.j = self.Next[self.i as usize];

                    while self.j != EMPTY {
                        // check if j and i have identical nonzero pattern

                        // check if i and j have the same Len and Elen

                        self.ok = if (self.Len[self.j as usize] == self.ln)
                            && (self.Elen[self.j as usize] == self.eln)
                        {
                            1
                        } else {
                            0
                        };
                        // skip the first element in the list (me)
                        self.p = self.Pe[self.j as usize] + 1;
                        while self.ok != 0 && self.p < self.Pe[self.j as usize] + self.ln {
                            if self.W[self.Iw[self.p as usize] as usize] != self.wflg {
                                self.ok = 0;
                            }
                            self.p += 1;
                        }
                        if self.ok != 0 {
                            // found it!  j can be absorbed into i

                            self.Pe[self.j as usize] = FLIP(self.i);
                            // both Nv [i] and Nv [j] are negated since they
                            // are in Lme, and the absolute values of each
                            // are the number of variables in i and j:
                            self.Nv[self.i as usize] += self.Nv[self.j as usize];
                            self.Nv[self.j as usize] = 0;
                            self.Elen[self.j as usize] = EMPTY;
                            // delete j from hash bucket

                            self.j = self.Next[self.j as usize];
                            self.Next[self.jlast as usize] = self.j;
                        } else {
                            // j cannot be absorbed into i
                            self.jlast = self.j;

                            self.j = self.Next[self.j as usize];
                        }
                    }

                    // no more variables can be absorbed into i
                    // go to next i in bucket and clear flag array

                    self.wflg += 1;
                    self.i = self.Next[self.i as usize];
                }
            }
            self.pme += 1;
        }

        // RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVARIABLES FROM ELEMENT (lines 1422-1473)

        self.p = self.pme1;
        self.nleft = self.n - self.nel;
        self.pme = self.pme1;
        while self.pme <= self.pme2 {
            self.i = self.Iw[self.pme as usize];

            self.nvi = -self.Nv[self.i as usize];

            if self.nvi > 0 {
                // i is a principal variable in Lme
                // restore Nv [i] to signify that i is principal
                self.Nv[self.i as usize] = self.nvi;

                // compute the external degree (add size of current element)

                self.deg = self.Degree[self.i as usize] + self.degme - self.nvi;
                self.deg = if self.deg < self.nleft - self.nvi {
                    self.deg
                } else {
                    self.nleft - self.nvi
                };

                // place the supervariable at the head of the degree list

                self.inext = self.Head[self.deg as usize];

                if self.inext != EMPTY {
                    self.Last[self.inext as usize] = self.i;
                }
                self.Next[self.i as usize] = self.inext;
                self.Last[self.i as usize] = EMPTY;
                self.Head[self.deg as usize] = self.i;

                // save the new degree, and find the minimum degree

                self.mindeg = if self.mindeg < self.deg {
                    self.mindeg
                } else {
                    self.deg
                };
                self.Degree[self.i as usize] = self.deg;

                // place the supervariable in the element pattern

                self.Iw[self.p as usize] = self.i;
                self.p += 1;
            }
            self.pme += 1;
        }

        // FINALIZE THE NEW ELEMENT (lines 1476-1524)

        self.Nv[self.me as usize] = self.nvpiv;
        // save the length of the list for the new element me
        self.Len[self.me as usize] = self.p - self.pme1;
        if self.Len[self.me as usize] == 0 {
            // there is nothing left of the current pivot element
            // it is a root of the assembly tree
            self.Pe[self.me as usize] = EMPTY;
            self.W[self.me as usize] = 0;
        }
        if self.elenme != 0 {
            // element was not constructed in place: deallocate part of
            // it since newly nonprincipal variables may have been removed
            self.pfree = self.p;
        }

        // The new element has nvpiv pivots and the size of the contribution
        // block for a multifrontal method is degme-by-degme, not including
        // the "dense" rows/columns.  If the "dense" rows/columns are included,
        // the frontal matrix is no larger than
        // (degme+ndense)-by-(degme+ndense).

        if self.Info.is_some() {
            self.f = self.nvpiv as f64;
            self.r = (self.degme + self.ndense) as f64;
            self.dmax = if self.dmax > self.f + self.r {
                self.dmax
            } else {
                self.f + self.r
            };

            // number of nonzeros in L (excluding the diagonal)
            self.lnzme = self.f * self.r + (self.f - 1.0) * self.f / 2.0;
            self.lnz += self.lnzme;

            // number of divide operations for LDL' and for LU
            self.ndiv += self.lnzme;

            // number of multiply-subtract pairs for LU
            self.s = self.f * self.r * self.r
                + self.r * (self.f - 1.0) * self.f
                + (self.f - 1.0) * self.f * (2.0 * self.f - 1.0) / 6.0;
            self.nms_lu += self.s;

            // number of multiply-subtract pairs for LDL'
            self.nms_ldl += (self.s + self.lnzme) / 2.0;
        }
    }
}
