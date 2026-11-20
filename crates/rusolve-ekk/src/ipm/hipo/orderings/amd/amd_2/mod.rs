//! HiGHS source: extern/amd/amd_2.c (lines 1-697, driver and initialisation)
//!
//! Translated from SuiteSparse AMD as vendored by HiGHS in `extern/amd` (BSD-3-Clause,
//! Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff).
//! See THIRD_PARTY_NOTICES.md.
//!
//! `amd_2`, the approximate minimum degree core. The C is one 1700-line
//! function; the Rust keeps every local in [`Amd2`] (same names, same types)
//! and splits the body at the C section banners into methods, called in the
//! C order: `init` (INITIALIZATIONS), then per pivot `construct` (GET PIVOT,
//! CONSTRUCT NEW ELEMENT), `degree` (COMPUTE W, DEGREE UPDATE), `superv`
//! (SUPERVARIABLE DETECTION, RESTORE DEGREE LISTS, FINALIZE NEW ELEMENT),
//! then `post` (DONE SELECTING PIVOTS, POST-ORDERING). Each method body is the
//! C section verbatim with `self.` locals; no reordering across sections.
//! `hash` is `amd_uint` so `hash % n` is the
//! unsigned remainder.

pub mod construct;
pub mod degree;
pub mod post;
pub mod superv;

use super::AMD_AGGRESSIVE;
use super::AMD_DEFAULT_AGGRESSIVE;
use super::AMD_DEFAULT_DENSE;
use super::AMD_DENSE;
use super::EMPTY;
use super::FLIP;
use super::amd_int;
use super::amd_int_max;
use super::amd_uint;

/// All arguments and locals of `Highs_amd_2` (amd_2.c lines 42-64 and the
/// declarations at lines 460-538).
pub struct Amd2<'a> {
    pub n: amd_int,
    pub Pe: &'a mut [amd_int],
    pub Iw: &'a mut [amd_int],
    pub Len: &'a mut [amd_int],
    pub iwlen: amd_int,
    pub pfree: amd_int,
    pub Nv: &'a mut [amd_int],
    pub Next: &'a mut [amd_int],
    pub Last: &'a mut [amd_int],
    pub Head: &'a mut [amd_int],
    pub Elen: &'a mut [amd_int],
    pub Degree: &'a mut [amd_int],
    pub W: &'a mut [amd_int],
    pub Control: Option<&'a [f64]>,
    pub Info: Option<&'a mut [f64]>,
    // LOCAL INTEGERS (lines 460-463)
    pub deg: amd_int,
    pub degme: amd_int,
    pub dext: amd_int,
    pub lemax: amd_int,
    pub e: amd_int,
    pub elenme: amd_int,
    pub eln: amd_int,
    pub i: amd_int,
    pub ilast: amd_int,
    pub inext: amd_int,
    pub j: amd_int,
    pub jlast: amd_int,
    pub jnext: amd_int,
    pub k: amd_int,
    pub knt1: amd_int,
    pub knt2: amd_int,
    pub knt3: amd_int,
    pub lenj: amd_int,
    pub ln: amd_int,
    pub me: amd_int,
    pub mindeg: amd_int,
    pub nel: amd_int,
    pub nleft: amd_int,
    pub nvi: amd_int,
    pub nvj: amd_int,
    pub nvpiv: amd_int,
    pub slenme: amd_int,
    pub wbig: amd_int,
    pub we: amd_int,
    pub wflg: amd_int,
    pub wnvi: amd_int,
    pub ok: amd_int,
    pub ndense: amd_int,
    pub ncmpa: amd_int,
    pub dense: amd_int,
    pub aggressive: amd_int,
    pub hash: amd_uint,
    // LOCAL DOUBLES (line 513)
    pub f: f64,
    pub r: f64,
    pub ndiv: f64,
    pub s: f64,
    pub nms_lu: f64,
    pub nms_ldl: f64,
    pub dmax: f64,
    pub alpha: f64,
    pub lnz: f64,
    pub lnzme: f64,
    // LOCAL "POINTERS" (line 538)
    pub p: amd_int,
    pub p1: amd_int,
    pub p2: amd_int,
    pub p3: amd_int,
    pub p4: amd_int,
    pub pdst: amd_int,
    pub pend: amd_int,
    pub pj: amd_int,
    pub pme: amd_int,
    pub pme1: amd_int,
    pub pme2: amd_int,
    pub pn: amd_int,
    pub psrc: amd_int,
}

/// HiGHS: amd_2.c::Highs_amd_2 (lines 42-64; algorithm notes 65-454; loop 690-697 and 1526-1530)
/// Builds [`Amd2`] with every local zero, then `init()`,
/// `while nel < n { construct(); degree(); superv(); }`, `post()`.
pub fn Highs_amd_2(
    n: amd_int,
    Pe: &mut [amd_int],
    Iw: &mut [amd_int],
    Len: &mut [amd_int],
    iwlen: amd_int,
    pfree: amd_int,
    Nv: &mut [amd_int],
    Next: &mut [amd_int],
    Last: &mut [amd_int],
    Head: &mut [amd_int],
    Elen: &mut [amd_int],
    Degree: &mut [amd_int],
    W: &mut [amd_int],
    Control: Option<&[f64]>,
    Info: Option<&mut [f64]>,
) {
    let mut a = Amd2 {
        n,
        Pe,
        Iw,
        Len,
        iwlen,
        pfree,
        Nv,
        Next,
        Last,
        Head,
        Elen,
        Degree,
        W,
        Control,
        Info,
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
        mindeg: 0,
        nel: 0,
        nleft: 0,
        nvi: 0,
        nvj: 0,
        nvpiv: 0,
        slenme: 0,
        wbig: 0,
        we: 0,
        wflg: 0,
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
    };
    a.init();
    while a.nel < a.n {
        a.construct();
        a.degree();
        a.superv();
    }
    a.post();
}

/// HiGHS: amd_2.c::Highs_amd_clear_flag (lines 22-35)
pub fn Highs_amd_clear_flag(
    wflg: amd_int,
    wbig: amd_int,
    W: &mut [amd_int],
    n: amd_int,
) -> amd_int {
    // Shared by the AMD_2 elimination steps.
    let mut wflg = wflg;
    if wflg < 2 || wflg >= wbig {
        for x in 0..n {
            if W[x as usize] != 0 {
                W[x as usize] = 1;
            }
        }
        wflg = 2;
    }
    // at this point, W [0..n-1] < wflg holds
    wflg
}

impl<'a> Amd2<'a> {
    /// HiGHS: amd_2.c::Highs_amd_2 INITIALIZATIONS (lines 556-689)
    pub fn init(&mut self) {
        // initialize output statistics
        self.lnz = 0.0;
        self.ndiv = 0.0;
        self.nms_lu = 0.0;
        self.nms_ldl = 0.0;
        self.dmax = 1.0;
        self.me = EMPTY;

        self.mindeg = 0;
        self.ncmpa = 0;
        self.nel = 0;
        self.lemax = 0;

        // get control parameters
        if let Some(control) = self.Control {
            self.alpha = control[AMD_DENSE];
            self.aggressive = if control[AMD_AGGRESSIVE] != 0.0 { 1 } else { 0 };
        } else {
            self.alpha = AMD_DEFAULT_DENSE;
            // C `#define AMD_DEFAULT_AGGRESSIVE 1` (an int); the skeleton
            // constant is f64 1.0, so cast back to the int local.
            self.aggressive = AMD_DEFAULT_AGGRESSIVE as amd_int;
        }
        // Note: if alpha is NaN, this is undefined:
        if self.alpha < 0.0 {
            // only remove completely dense rows/columns
            self.dense = self.n - 2;
        } else {
            self.dense = (self.alpha * (self.n as f64).sqrt()) as amd_int;
        }
        // PORT NOTE: AMD MAX is a C ternary, not libm fmax
        self.dense = if 16 > self.dense { 16 } else { self.dense };
        // PORT NOTE: AMD MIN is a C ternary, not libm fmin
        self.dense = if self.n < self.dense {
            self.n
        } else {
            self.dense
        };

        self.i = 0;
        while self.i < self.n {
            let i = self.i;
            self.Last[i as usize] = EMPTY;
            self.Head[i as usize] = EMPTY;
            self.Next[i as usize] = EMPTY;
            // if separate Hhead array is used for hash buckets: *
            // Hhead [i] = EMPTY ;
            self.Nv[i as usize] = 1;
            self.W[i as usize] = 1;
            self.Elen[i as usize] = 0;
            self.Degree[i as usize] = self.Len[i as usize];
            self.i += 1;
        }

        // initialize wflg
        self.wbig = amd_int_max - self.n;
        self.wflg = Highs_amd_clear_flag(0, self.wbig, self.W, self.n);

        // ---------------------------------------------------------------------
        // initialize degree lists and eliminate dense and empty rows
        // ---------------------------------------------------------------------

        self.ndense = 0;

        self.i = 0;
        while self.i < self.n {
            let i = self.i;
            self.deg = self.Degree[i as usize];

            if self.deg == 0 {
                // -------------------------------------------------------------
                // we have a variable that can be eliminated at once because
                // there is no off-diagonal non-zero in its row.  Note that
                // Nv [i] = 1 for an empty variable i.  It is treated just
                // the same as an eliminated element i.
                // -------------------------------------------------------------

                self.Elen[i as usize] = FLIP(1);
                self.nel += 1;
                self.Pe[i as usize] = EMPTY;
                self.W[i as usize] = 0;
            } else if self.deg > self.dense {
                // -------------------------------------------------------------
                // Dense variables are not treated as elements, but as unordered,
                // non-principal variables that have no parent.  They do not take
                // part in the postorder, since Nv [i] = 0.  Note that the Fortran
                // version does not have this option.
                // -------------------------------------------------------------

                self.ndense += 1;
                self.Nv[i as usize] = 0; // do not postorder this node
                self.Elen[i as usize] = EMPTY;
                self.nel += 1;
                self.Pe[i as usize] = EMPTY;
            } else {
                // -------------------------------------------------------------
                // place i in the degree list corresponding to its degree
                // -------------------------------------------------------------

                self.inext = self.Head[self.deg as usize];

                if self.inext != EMPTY {
                    self.Last[self.inext as usize] = i;
                }
                self.Next[i as usize] = self.inext;
                self.Head[self.deg as usize] = i;
            }
            self.i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn clear_flag_resets_when_big() {
        let mut w = [5, 0, 7];
        let wflg = Highs_amd_clear_flag(10, 10, &mut w, 3);
        assert_eq!(wflg, 2);
        assert_eq!(w, [1, 0, 1]);
    }
}
