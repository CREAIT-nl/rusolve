//! HiGHS source: extern/metis/libmetis/wspace.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Workspace (`gk_mcore_t`) replacement. The C hands out uninitialised slices
//! of one core buffer and frees them at `WCOREPOP`; correct METIS code writes
//! every workspace entry before reading it, so a fresh zeroed `Vec` per call is
//! observably identical.

use super::{ctrl_t, graph_t, idx_t, ikv_t, real_t};

/// HiGHS: wspace.c::AllocateWorkSpace (lines 17-28)
pub fn AllocateWorkSpace(ctrl: &mut ctrl_t, graph: &graph_t) {
    ctrl.mcore = 0;
}

/// HiGHS: wspace.c::FreeWorkSpace (lines 33-36)
pub fn FreeWorkSpace(ctrl: &mut ctrl_t) {
    ctrl.mcore = 0;
}

/// HiGHS: wspace.c::wspacepush (lines 52-55), `WCOREPUSH`
pub fn wspacepush(ctrl: &mut ctrl_t) {
    ctrl.mcore += 1;
}

/// HiGHS: wspace.c::wspacepop (lines 61-64), `WCOREPOP`
pub fn wspacepop(ctrl: &mut ctrl_t) {
    assert!(ctrl.mcore > 0, "WCOREPOP without WCOREPUSH");
    ctrl.mcore -= 1;
}

/// HiGHS: wspace.c::iwspacemalloc (lines 70-73)
pub fn iwspacemalloc(ctrl: &mut ctrl_t, n: idx_t) -> Vec<idx_t> {
    vec![0; n.max(0) as usize]
}

/// HiGHS: wspace.c::rwspacemalloc (lines 79-82)
pub fn rwspacemalloc(ctrl: &mut ctrl_t, n: idx_t) -> Vec<real_t> {
    vec![0.0; n.max(0) as usize]
}

/// HiGHS: wspace.c::ikvwspacemalloc (lines 88-91)
pub fn ikvwspacemalloc(ctrl: &mut ctrl_t, n: idx_t) -> Vec<ikv_t> {
    vec![ikv_t::default(); n.max(0) as usize]
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn push_pop_depth() {
        let mut ctrl = ctrl_t::default();
        wspacepush(&mut ctrl);
        let v = iwspacemalloc(&mut ctrl, 3);
        assert_eq!(v, vec![0, 0, 0]);
        wspacepop(&mut ctrl);
        assert_eq!(ctrl.mcore, 0);
    }
}
