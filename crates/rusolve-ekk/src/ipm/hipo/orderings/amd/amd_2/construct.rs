//! HiGHS source: extern/amd/amd_2.c (lines 698-997)
//!
//! Translated from SuiteSparse AMD as vendored by HiGHS in `extern/amd` (BSD-3-Clause,
//! Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff).
//! See THIRD_PARTY_NOTICES.md.
//!
//! GET PIVOT OF MINIMUM DEGREE and CONSTRUCT NEW ELEMENT (including the
//! in-place garbage collection of `Iw`, which counts `ncmpa`).

use super::super::{EMPTY, FLIP};
use super::Amd2;

impl<'a> Amd2<'a> {
    /// HiGHS: amd_2.c::Highs_amd_2 GET PIVOT + CONSTRUCT NEW ELEMENT (lines 698-997)
    pub fn construct(&mut self) {
        // -----------------------------------------------------------------
        // find next supervariable for elimination
        // -----------------------------------------------------------------

        self.deg = self.mindeg;
        while self.deg < self.n {
            self.me = self.Head[self.deg as usize];
            if self.me != EMPTY {
                break;
            }
            self.deg += 1;
        }
        self.mindeg = self.deg;

        // -----------------------------------------------------------------
        // remove chosen variable from link list
        // -----------------------------------------------------------------

        self.inext = self.Next[self.me as usize];

        if self.inext != EMPTY {
            self.Last[self.inext as usize] = EMPTY;
        }
        self.Head[self.deg as usize] = self.inext;

        // -----------------------------------------------------------------
        // me represents the elimination of pivots nel to nel+Nv[me]-1.
        // place me itself as the first in this set.
        // -----------------------------------------------------------------

        self.elenme = self.Elen[self.me as usize];
        self.nvpiv = self.Nv[self.me as usize];

        self.nel += self.nvpiv;

        // =========================================================================
        // CONSTRUCT NEW ELEMENT
        // =========================================================================

        // -----------------------------------------------------------------
        // At this point, me is the pivotal supervariable.  It will be
        // converted into the current element.  Scan list of the pivotal
        // supervariable, me, setting tree pointers and constructing new list
        // of supervariables for the new element, me.  p is a pointer to the
        // current position in the old list.
        // -----------------------------------------------------------------

        // flag the variable "me" as being in Lme by negating Nv [me]
        self.Nv[self.me as usize] = -self.nvpiv;
        self.degme = 0;

        if self.elenme == 0 {
            // -------------------------------------------------------------
            // construct the new element in place
            // -------------------------------------------------------------

            self.pme1 = self.Pe[self.me as usize];
            self.pme2 = self.pme1 - 1;

            self.p = self.pme1;
            while self.p < self.pme1 + self.Len[self.me as usize] {
                self.i = self.Iw[self.p as usize];

                self.nvi = self.Nv[self.i as usize];
                if self.nvi > 0 {
                    // -----------------------------------------------------
                    // i is a principal variable not yet placed in Lme.
                    // store i in new list
                    // -----------------------------------------------------

                    // flag i as being in Lme by negating Nv [i]
                    self.degme += self.nvi;
                    self.Nv[self.i as usize] = -self.nvi;
                    self.pme2 += 1;
                    self.Iw[self.pme2 as usize] = self.i;

                    // -----------------------------------------------------
                    // remove variable i from degree list.
                    // -----------------------------------------------------

                    self.ilast = self.Last[self.i as usize];
                    self.inext = self.Next[self.i as usize];

                    if self.inext != EMPTY {
                        self.Last[self.inext as usize] = self.ilast;
                    }
                    if self.ilast != EMPTY {
                        self.Next[self.ilast as usize] = self.inext;
                    } else {
                        // i is at the head of the degree list

                        self.Head[self.Degree[self.i as usize] as usize] = self.inext;
                    }
                }
                self.p += 1;
            }
        } else {
            // -------------------------------------------------------------
            // construct the new element in empty space, Iw [pfree ...]
            // -------------------------------------------------------------

            self.p = self.Pe[self.me as usize];
            self.pme1 = self.pfree;
            self.slenme = self.Len[self.me as usize] - self.elenme;

            self.knt1 = 1;
            while self.knt1 <= self.elenme + 1 {
                if self.knt1 > self.elenme {
                    // search the supervariables in me.
                    self.e = self.me;
                    self.pj = self.p;
                    self.ln = self.slenme;
                } else {
                    // search the elements in me.
                    self.e = self.Iw[self.p as usize];
                    self.p += 1;

                    self.pj = self.Pe[self.e as usize];
                    self.ln = self.Len[self.e as usize];
                }

                // ---------------------------------------------------------
                // search for different supervariables and add them to the
                // new list, compressing when necessary. this loop is
                // executed once for each element in the list and once for
                // all the supervariables in the list.
                // ---------------------------------------------------------

                self.knt2 = 1;
                while self.knt2 <= self.ln {
                    self.i = self.Iw[self.pj as usize];
                    self.pj += 1;

                    self.nvi = self.Nv[self.i as usize];

                    if self.nvi > 0 {
                        // -------------------------------------------------
                        // compress Iw, if necessary
                        // -------------------------------------------------

                        if self.pfree >= self.iwlen {
                            // prepare for compressing Iw by adjusting pointers
                            // and lengths so that the lists being searched in
                            // the inner and outer loops contain only the
                            // remaining entries.

                            self.Pe[self.me as usize] = self.p;
                            self.Len[self.me as usize] -= self.knt1;
                            // check if nothing left of supervariable me
                            if self.Len[self.me as usize] == 0 {
                                self.Pe[self.me as usize] = EMPTY;
                            }
                            self.Pe[self.e as usize] = self.pj;
                            self.Len[self.e as usize] = self.ln - self.knt2;
                            // nothing left of element e
                            if self.Len[self.e as usize] == 0 {
                                self.Pe[self.e as usize] = EMPTY;
                            }

                            self.ncmpa += 1; // one more garbage collection

                            // store first entry of each object in Pe
                            // FLIP the first entry in each object
                            self.j = 0;
                            while self.j < self.n {
                                self.pn = self.Pe[self.j as usize];
                                if self.pn >= 0 {
                                    self.Pe[self.j as usize] = self.Iw[self.pn as usize];
                                    self.Iw[self.pn as usize] = FLIP(self.j);
                                }
                                self.j += 1;
                            }

                            // psrc/pdst point to source/destination
                            self.psrc = 0;
                            self.pdst = 0;
                            self.pend = self.pme1 - 1;

                            while self.psrc <= self.pend {
                                // search for next FLIP'd entry
                                self.j = FLIP(self.Iw[self.psrc as usize]);
                                self.psrc += 1;
                                if self.j >= 0 {
                                    self.Iw[self.pdst as usize] = self.Pe[self.j as usize];
                                    self.Pe[self.j as usize] = self.pdst;
                                    self.pdst += 1;
                                    self.lenj = self.Len[self.j as usize];
                                    // copy from source to destination
                                    self.knt3 = 0;
                                    while self.knt3 <= self.lenj - 2 {
                                        self.Iw[self.pdst as usize] = self.Iw[self.psrc as usize];
                                        self.pdst += 1;
                                        self.psrc += 1;
                                        self.knt3 += 1;
                                    }
                                }
                            }

                            // move the new partially-constructed element
                            self.p1 = self.pdst;
                            self.psrc = self.pme1;
                            while self.psrc < self.pfree {
                                self.Iw[self.pdst as usize] = self.Iw[self.psrc as usize];
                                self.pdst += 1;
                                self.psrc += 1;
                            }
                            self.pme1 = self.p1;
                            self.pfree = self.pdst;
                            self.pj = self.Pe[self.e as usize];
                            self.p = self.Pe[self.me as usize];
                        }

                        // -------------------------------------------------
                        // i is a principal variable not yet placed in Lme
                        // store i in new list
                        // -------------------------------------------------

                        // flag i as being in Lme by negating Nv [i]
                        self.degme += self.nvi;
                        self.Nv[self.i as usize] = -self.nvi;
                        self.Iw[self.pfree as usize] = self.i;
                        self.pfree += 1;

                        // -------------------------------------------------
                        // remove variable i from degree link list
                        // -------------------------------------------------

                        self.ilast = self.Last[self.i as usize];
                        self.inext = self.Next[self.i as usize];

                        if self.inext != EMPTY {
                            self.Last[self.inext as usize] = self.ilast;
                        }
                        if self.ilast != EMPTY {
                            self.Next[self.ilast as usize] = self.inext;
                        } else {
                            // i is at the head of the degree list

                            self.Head[self.Degree[self.i as usize] as usize] = self.inext;
                        }
                    }
                    self.knt2 += 1;
                }

                if self.e != self.me {
                    // set tree pointer and flag to indicate element e is
                    // absorbed into new element me (the parent of e is me)

                    self.Pe[self.e as usize] = FLIP(self.me);
                    self.W[self.e as usize] = 0;
                }
                self.knt1 += 1;
            }

            self.pme2 = self.pfree - 1;
        }

        // -----------------------------------------------------------------
        // me has now been converted into an element in Iw [pme1..pme2]
        // -----------------------------------------------------------------

        // degme holds the external degree of new element
        self.Degree[self.me as usize] = self.degme;
        self.Pe[self.me as usize] = self.pme1;
        self.Len[self.me as usize] = self.pme2 - self.pme1 + 1;

        self.Elen[self.me as usize] = FLIP(self.nvpiv + self.degme);
        // FLIP (Elen (me)) is now the degree of pivot (including
        // diagonal part).

        // -----------------------------------------------------------------
        // make sure that wflg is not too large.
        // -----------------------------------------------------------------

        // With the current value of wflg, wflg+n must not cause integer
        // overflow

        self.wflg = super::Highs_amd_clear_flag(self.wflg, self.wbig, self.W, self.n);
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests, hand-derived from amd_2.c 698-997.
    use super::*;
    use crate::ipm::hipo::orderings::amd::amd_int;

    #[allow(clippy::too_many_arguments)]
    fn state<'a>(
        n: amd_int,
        pe: &'a mut [amd_int],
        iw: &'a mut [amd_int],
        len: &'a mut [amd_int],
        nv: &'a mut [amd_int],
        next: &'a mut [amd_int],
        last: &'a mut [amd_int],
        head: &'a mut [amd_int],
        elen: &'a mut [amd_int],
        degree: &'a mut [amd_int],
        w: &'a mut [amd_int],
        pfree: amd_int,
        iwlen: amd_int,
        mindeg: amd_int,
        wflg: amd_int,
        wbig: amd_int,
    ) -> Amd2<'a> {
        Amd2 {
            n,
            Pe: pe,
            Iw: iw,
            Len: len,
            iwlen,
            pfree,
            Nv: nv,
            Next: next,
            Last: last,
            Head: head,
            Elen: elen,
            Degree: degree,
            W: w,
            Control: None,
            Info: None,
            deg: 0,
            degme: 0,
            dext: 0,
            lemax: 0,
            e: 0,
            elenme: 0,
            eln: 0,
            i: 0,
            ilast: 0,
            inext: 0,
            j: 0,
            jlast: 0,
            jnext: 0,
            k: 0,
            knt1: 0,
            knt2: 0,
            knt3: 0,
            lenj: 0,
            ln: 0,
            me: 0,
            mindeg,
            nel: 0,
            nleft: 0,
            nvi: 0,
            nvj: 0,
            nvpiv: 0,
            slenme: 0,
            wbig,
            we: 0,
            wflg,
            wnvi: 0,
            ok: 0,
            ndense: 0,
            ncmpa: 0,
            dense: 0,
            aggressive: 0,
            hash: 0,
            f: 0.0,
            r: 0.0,
            ndiv: 0.0,
            s: 0.0,
            nms_lu: 0.0,
            nms_ldl: 0.0,
            dmax: 0.0,
            alpha: 0.0,
            lnz: 0.0,
            lnzme: 0.0,
            p: 0,
            p1: 0,
            p2: 0,
            p3: 0,
            p4: 0,
            pdst: 0,
            pend: 0,
            pj: 0,
            pme: 0,
            pme1: 0,
            pme2: 0,
            pn: 0,
            psrc: 0,
        }
    }

    /// Pivot on a leaf supervariable (elenme == 0, in-place branch). Variables 0..3,
    /// me = 1 at degree 1 with adjacency {2}; var 2 is also at degree 1 and is absorbed.
    #[test]
    fn construct_leaf_pivot_in_place() {
        let mut pe = vec![EMPTY, 0, EMPTY];
        let mut iw = vec![2, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let mut len = vec![0, 1, 0];
        let mut nv = vec![1, 1, 1];
        let mut next = vec![EMPTY, 2, EMPTY];
        let mut last = vec![EMPTY, EMPTY, 1];
        let mut head = vec![EMPTY, 1, EMPTY];
        let mut elen = vec![0, 0, 0];
        let mut degree = vec![0, 1, 1];
        let mut w = vec![0, 0, 0];
        let mut s = state(
            3,
            &mut pe,
            &mut iw,
            &mut len,
            &mut nv,
            &mut next,
            &mut last,
            &mut head,
            &mut elen,
            &mut degree,
            &mut w,
            3,
            10,
            1,
            2,
            100,
        );
        s.construct();
        assert_eq!(s.me, 1);
        assert_eq!(s.nvpiv, 1);
        assert_eq!(s.nel, 1);
        assert_eq!(s.mindeg, 1);
        assert_eq!(s.degme, 1);
        assert_eq!(s.wflg, 2);
        assert_eq!(s.Nv, [1, -1, -1]);
        assert_eq!(s.Iw[0], 2);
        assert_eq!(s.Pe[1], 0);
        assert_eq!(s.Len[1], 1);
        assert_eq!(s.Degree[1], 1);
        assert_eq!(s.Elen[1], FLIP(2));
        assert_eq!(s.Head[1], EMPTY);
    }

    /// elenme > 0 branch with an Iw garbage collection (pfree == iwlen). me = 0
    /// is adjacent to element 2 (list {0, 1}) and variable 1; var 1 is adjacent to
    /// element 2 and var 0. Traced by hand through amd_2.c 815-958.
    #[test]
    fn construct_from_element_with_compression() {
        // Iw: element 2 at 0..2 = [0, 1]; var 0 at 2..4 = [2, 1]; var 1 at 4..6 = [2, 0]
        let mut pe = vec![2, 4, 0];
        let mut iw = vec![0, 1, 2, 1, 2, 0];
        let mut len = vec![2, 2, 2];
        let mut nv = vec![1, 1, 1];
        let mut next = vec![1, EMPTY, EMPTY];
        let mut last = vec![EMPTY, 0, EMPTY];
        let mut head = vec![EMPTY, 0, EMPTY];
        let mut elen = vec![1, 1, EMPTY];
        let mut degree = vec![1, 1, 0];
        let mut w = vec![1, 1, 5];
        let mut s = state(
            3,
            &mut pe,
            &mut iw,
            &mut len,
            &mut nv,
            &mut next,
            &mut last,
            &mut head,
            &mut elen,
            &mut degree,
            &mut w,
            6,
            6,
            0,
            2,
            100,
        );
        s.construct();
        assert_eq!(s.me, 0);
        assert_eq!(s.mindeg, 1);
        assert_eq!(s.nel, 1);
        assert_eq!(s.ncmpa, 1);
        assert_eq!(s.degme, 1);
        assert_eq!(s.pme1, 3);
        assert_eq!(s.pme2, 3);
        assert_eq!(s.pfree, 4);
        assert_eq!(s.wflg, 2);
        assert_eq!(s.Iw, [1, 2, 0, 1, FLIP(1), 0]);
        assert_eq!(s.Pe, [3, 1, FLIP(0)]);
        assert_eq!(s.Len, [1, 2, 0]);
        assert_eq!(s.Nv, [-1, -1, 1]);
        assert_eq!(s.Elen, [FLIP(2), 1, EMPTY]);
        assert_eq!(s.Degree, [1, 1, 0]);
        assert_eq!(s.W, [1, 1, 0]);
        assert_eq!(s.Head, [EMPTY, EMPTY, EMPTY]);
    }
}
