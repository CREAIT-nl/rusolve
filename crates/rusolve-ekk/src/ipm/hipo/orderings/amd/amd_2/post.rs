//! HiGHS source: extern/amd/amd_2.c (lines 1527-1761)
//!
//! Translated from SuiteSparse AMD as vendored by HiGHS in `extern/amd` (BSD-3-Clause,
//! Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff).
//! See THIRD_PARTY_NOTICES.md.
//!
//! DONE SELECTING PIVOTS (dense-row statistics into `Info`) and POST-ORDERING
//! (path compression of the absorbed variables, `amd_postorder` of the
//! element tree, output permutation in `Last` / `Next`).

use super::super::AMD_DMAX;
use super::super::AMD_LNZ;
use super::super::AMD_NCMPA;
use super::super::AMD_NDENSE;
use super::super::AMD_NDIV;
use super::super::AMD_NMULTSUBS_LDL;
use super::super::AMD_NMULTSUBS_LU;
use super::super::AMD_OK;
use super::super::AMD_STATUS;
use super::super::EMPTY;
use super::super::FLIP;
use super::super::Highs_amd_postorder;
use super::super::amd_int;
use super::Amd2;

impl<'a> Amd2<'a> {
    /// HiGHS: amd_2.c::Highs_amd_2 DONE SELECTING PIVOTS + POST-ORDERING (lines 1527-1761)
    pub fn post(&mut self) {
        if let Some(info) = self.Info.as_deref_mut() {
            // count the work to factorize the ndense-by-ndense submatrix
            self.f = self.ndense as f64;
            // PORT NOTE: AMD MAX is a C ternary, not libm fmax
            self.dmax = if self.dmax > (self.ndense as f64) {
                self.dmax
            } else {
                self.ndense as f64
            };

            // number of nonzeros in L (excluding the diagonal)
            self.lnzme = (self.f - 1.0) * self.f / 2.0;
            self.lnz += self.lnzme;

            // number of divide operations for LDL' and for LU
            self.ndiv += self.lnzme;

            // number of multiply-subtract pairs for LU
            self.s = (self.f - 1.0) * self.f * (2.0 * self.f - 1.0) / 6.0;
            self.nms_lu += self.s;

            // number of multiply-subtract pairs for LDL'
            self.nms_ldl += (self.s + self.lnzme) / 2.0;

            // number of nz's in L (excl. diagonal)
            info[AMD_LNZ] = self.lnz;

            // number of divide ops for LU and LDL'
            info[AMD_NDIV] = self.ndiv;

            // number of multiply-subtract pairs for LDL'
            info[AMD_NMULTSUBS_LDL] = self.nms_ldl;

            // number of multiply-subtract pairs for LU
            info[AMD_NMULTSUBS_LU] = self.nms_lu;

            // number of "dense" rows/columns
            info[AMD_NDENSE] = self.ndense as f64;

            // largest front is dmax-by-dmax
            info[AMD_DMAX] = self.dmax;

            // number of garbage collections in AMD
            info[AMD_NCMPA] = self.ncmpa as f64;

            // successful ordering
            info[AMD_STATUS] = AMD_OK as f64;
        }

        // restore Pe
        self.i = 0;
        while self.i < self.n {
            self.Pe[self.i as usize] = FLIP(self.Pe[self.i as usize]);
            self.i += 1;
        }

        // restore Elen, for output information, and for postordering
        self.i = 0;
        while self.i < self.n {
            self.Elen[self.i as usize] = FLIP(self.Elen[self.i as usize]);
            self.i += 1;
        }

        // Now the parent of j is Pe [j], or EMPTY if j is a root.  Elen [e] > 0
        // is the size of element e.  Elen [i] is EMPTY for unordered variable i.

        // compress the paths of the variables
        self.i = 0;
        while self.i < self.n {
            let i = self.i;
            if self.Nv[i as usize] == 0 {
                // -------------------------------------------------------------
                // i is an un-ordered row.  Traverse the tree from i until
                // reaching an element, e.  The element, e, was the principal
                // supervariable of i and all nodes in the path from i to when e
                // was selected as pivot.
                // -------------------------------------------------------------

                self.j = self.Pe[i as usize];

                if self.j == EMPTY {
                    // Skip a dense variable.  It has no parent.
                    // (the C for-loop increment runs on continue)
                    self.i += 1;
                    continue;
                }

                // while (j is a variable)
                while self.Nv[self.j as usize] == 0 {
                    self.j = self.Pe[self.j as usize];
                }
                // got to an element e
                self.e = self.j;

                // -------------------------------------------------------------
                // traverse the path again from i to e, and compress the path
                // (all nodes point to e).  Path compression allows this code to
                // compute in O(n) time.
                // -------------------------------------------------------------

                self.j = i;
                // while (j is a variable)
                while self.Nv[self.j as usize] == 0 {
                    self.jnext = self.Pe[self.j as usize];

                    self.Pe[self.j as usize] = self.e;
                    self.j = self.jnext;
                }
            }
            self.i += 1;
        }

        // postorder the assembly tree
        Highs_amd_postorder(
            self.n, self.Pe, self.Nv, self.Elen, self.W, // output order
            self.Head, self.Next, self.Last,
        );

        // compute output permutation and inverse permutation

        // W [e] = k means that element e is the kth element in the new
        // order.  e is in the range 0 to n-1, and k is in the range 0 to
        // the number of elements.  Use Head for inverse order.

        self.k = 0;
        while self.k < self.n {
            self.Head[self.k as usize] = EMPTY;
            self.Next[self.k as usize] = EMPTY;
            self.k += 1;
        }
        self.e = 0;
        while self.e < self.n {
            self.k = self.W[self.e as usize];

            if self.k != EMPTY {
                self.Head[self.k as usize] = self.e;
            }
            self.e += 1;
        }

        // construct output inverse permutation in Next,
        // and permutation in Last
        self.nel = 0;
        self.k = 0;
        while self.k < self.n {
            self.e = self.Head[self.k as usize];
            if self.e == EMPTY {
                break;
            }

            self.Next[self.e as usize] = self.nel;
            self.nel += self.Nv[self.e as usize];
            self.k += 1;
        }

        // order non-principal variables (dense, & those merged into supervar's)
        self.i = 0;
        while self.i < self.n {
            let i = self.i;
            if self.Nv[i as usize] == 0 {
                self.e = self.Pe[i as usize];

                if self.e != EMPTY {
                    // This is an unordered variable that was merged
                    // into element e via supernode detection or mass
                    // elimination of i when e became the pivot element.
                    // Place i in order just before e.

                    self.Next[i as usize] = self.Next[self.e as usize];
                    self.Next[self.e as usize] += 1;
                } else {
                    // This is a dense unordered variable, with no parent.
                    // Place it last in the output order.
                    self.Next[i as usize] = self.nel;
                    self.nel += 1;
                }
            }
            self.i += 1;
        }

        self.i = 0;
        while self.i < self.n {
            self.k = self.Next[self.i as usize];

            self.Last[self.k as usize] = self.i;
            self.i += 1;
        }
    }
}
