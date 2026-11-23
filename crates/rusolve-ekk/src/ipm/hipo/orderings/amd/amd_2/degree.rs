//! HiGHS source: extern/amd/amd_2.c (lines 998-1280)
//!
//! Translated from SuiteSparse AMD as vendored by HiGHS in `extern/amd` (BSD-3-Clause,
//! Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff).
//! See THIRD_PARTY_NOTICES.md.
//!
//! COMPUTE (W [e] - wflg) = |Le\Lme| FOR ALL ELEMENTS and DEGREE UPDATE AND
//! ELEMENT ABSORPTION (aggressive absorption when `aggressive != 0`).

use super::super::{EMPTY, FLIP, amd_int, amd_uint};
use super::Amd2;
use super::Highs_amd_clear_flag;

impl<'a> Amd2<'a> {
    /// HiGHS: amd_2.c::Highs_amd_2 COMPUTE W + DEGREE UPDATE (lines 998-1280)
    ///
    /// PORT NOTE (rusolve): the same statements and
    /// order as the C, but the scratch scalars live in locals instead of
    /// `self` fields and are written back at the end. With every `p`, `e`,
    /// `we` a field behind `&mut self` the compiler reloads and stores them
    /// through memory on each iteration; on rail4284 this loop was 3640 of
    /// 3670 profile samples and the port's analyse ran 1.9x slower than the
    /// C. The result is bit-identical (same integer operations).
    pub fn degree(&mut self) {
        let n = self.n;
        let me = self.me;
        let pme1 = self.pme1;
        let pme2 = self.pme2;
        let wflg = self.wflg;
        let aggressive = self.aggressive;
        let mut degme = self.degme;
        let mut nvpiv = self.nvpiv;
        let mut nel = self.nel;
        let iw = &mut *self.Iw;
        let pe = &mut *self.Pe;
        let elen = &mut *self.Elen;
        let nv = &mut *self.Nv;
        let w = &mut *self.W;
        let degree = &mut *self.Degree;
        let len = &mut *self.Len;
        let head = &mut *self.Head;
        let next = &mut *self.Next;
        let last = &mut *self.Last;

        // COMPUTE (W [e] - wflg) = |Le\Lme| FOR ALL ELEMENTS (lines 998-1051)
        // Scan 1: compute the external degrees of previous elements with
        // respect to the current element.
        for pme in pme1..=pme2 {
            let i = iw[pme as usize];
            let eln = elen[i as usize];
            if eln > 0 {
                // note that Nv [i] has been negated to denote i in Lme:
                let nvi = -nv[i as usize];
                let wnvi = wflg - nvi;
                let p0 = pe[i as usize];
                for p in p0..p0 + eln {
                    let e = iw[p as usize];
                    let mut we = w[e as usize];
                    if we >= wflg {
                        // unabsorbed element e has been seen in this loop
                        we -= nvi;
                    } else if we != 0 {
                        // e is an unabsorbed element
                        // this is the first we have seen e in all of Scan 1
                        we = degree[e as usize] + wnvi;
                    }
                    w[e as usize] = we;
                }
            }
        }

        // DEGREE UPDATE AND ELEMENT ABSORPTION (lines 1054-1266)
        // Scan 2: for each i in Lme, sum up the degree of Lme (which is
        // degme), plus the sum of the external degrees of each Le for the
        // elements e appearing within i, plus the supervariables in i.
        // Place i in hash list.
        for pme in pme1..=pme2 {
            let i = iw[pme as usize];
            let p1 = pe[i as usize];
            let p2 = p1 + elen[i as usize] - 1;
            let mut pn = p1;
            // hash is unsigned Int in C: `hash += e` wraps, hence wrapping_add.
            let mut hash: amd_uint = 0;
            let mut deg: amd_int = 0;

            // scan the element list associated with supervariable i (lines 1077-1131)
            // UMFPACK/MA38-style approximate degree:
            if aggressive != 0 {
                for p in p1..=p2 {
                    let e = iw[p as usize];
                    let we = w[e as usize];
                    if we != 0 {
                        // e is an unabsorbed element
                        // dext = | Le \ Lme |
                        let dext = we - wflg;
                        if dext > 0 {
                            deg += dext;
                            iw[pn as usize] = e;
                            pn += 1;
                            hash = hash.wrapping_add(e as amd_uint);
                        } else {
                            // external degree of e is zero, absorb e into me
                            pe[e as usize] = FLIP(me);
                            w[e as usize] = 0;
                        }
                    }
                }
            } else {
                for p in p1..=p2 {
                    let e = iw[p as usize];
                    let we = w[e as usize];
                    if we != 0 {
                        // e is an unabsorbed element
                        let dext = we - wflg;
                        deg += dext;
                        iw[pn as usize] = e;
                        pn += 1;
                        hash = hash.wrapping_add(e as amd_uint);
                    }
                }
            }

            // count the number of elements in i (including me):
            elen[i as usize] = pn - p1 + 1;

            // scan the supervariables in the list associated with i (lines 1136-1160)
            let p3 = pn;
            let p4 = p1 + len[i as usize];
            let mut p = p2 + 1;
            while p < p4 {
                let j = iw[p as usize];
                let nvj = nv[j as usize];
                if nvj > 0 {
                    // j is unabsorbed, and not in Lme.
                    // add to degree and add to new list
                    deg += nvj;
                    iw[pn as usize] = j;
                    pn += 1;
                    hash = hash.wrapping_add(j as amd_uint);
                }
                p += 1;
            }

            // update the degree and check for mass elimination (lines 1162-1265)
            // with aggressive absorption, deg==0 is identical to the
            // Elen [i] == 1 && p3 == pn test, below.
            if elen[i as usize] == 1 && p3 == pn {
                // mass elimination (lines 1173-1206)
                pe[i as usize] = FLIP(me);
                let nvi = -nv[i as usize];
                degme -= nvi;
                nvpiv += nvi;
                nel += nvi;
                nv[i as usize] = 0;
                elen[i as usize] = EMPTY;
            } else {
                // update the upper-bound degree of i (lines 1210-1265)
                // the following degree does not yet include the size
                // of the current element, which is added later:
                // MIN is a C ternary: if a < b { a } else { b }
                degree[i as usize] = if degree[i as usize] < deg {
                    degree[i as usize]
                } else {
                    deg
                };

                // add me to the list for i
                // move first supervariable to end of list
                iw[pn as usize] = iw[p3 as usize];
                // move first element to end of element part of list
                iw[p3 as usize] = iw[p1 as usize];
                // add new element, me, to front of list.
                iw[p1 as usize] = me;
                // store the new length of the list in Len [i]
                len[i as usize] = pn - p1 + 1;

                // place in hash bucket.  Save hash key of i in Last [i].
                // NOTE: this can fail if hash is negative, because the ANSI C
                // standard does not define a % b when a and/or b are negative.
                // That's why hash is defined as an unsigned Int, to avoid this
                // problem.
                hash %= n as amd_uint;

                // if the Hhead array is not used:
                let j = head[hash as usize];
                if j <= EMPTY {
                    // degree list is empty, hash head is FLIP (j)
                    next[i as usize] = FLIP(j);
                    head[hash as usize] = FLIP(i);
                } else {
                    // degree list is not empty, use Last [Head [hash]] as
                    // hash head.
                    next[i as usize] = last[j as usize];
                    last[j as usize] = i;
                }
                last[i as usize] = hash as amd_int;
            }
        }

        degree[me as usize] = degme;
        self.degme = degme;
        self.nvpiv = nvpiv;
        self.nel = nel;

        // Clear the counter array, W [...], by incrementing wflg. (lines 1270-1278)
        // make sure that wflg+n does not cause integer overflow
        // MAX is a C ternary: if a > b { a } else { b }
        self.lemax = if self.lemax > degme {
            self.lemax
        } else {
            degme
        };
        self.wflg += self.lemax;
        self.wflg = Highs_amd_clear_flag(self.wflg, self.wbig, self.W, self.n);
        //  at this point, W [0..n-1] < wflg holds
    }
}
