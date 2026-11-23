//! HiGHS source: extern/metis/metis.h, extern/metis/libmetis/{defs.h,struct.h,macros.h,gklib_defs.h}, extern/metis/GKlib/{gk_macros.h,gk_struct.h}
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! METIS 5.2.1 nested-dissection ordering, trimmed to the call graph of
//! `Highs_METIS_NodeND` (the only METIS entry point HiPO uses, through
//! `HighsExtras::metis::nodeND`). Only `METIS_OP_OMETIS` exists in the vendored
//! copy. Type widths follow `metis.h` as configured by HiGHS: `idx_t` is
//! [`crate::HighsInt`] (`IDXTYPEWIDTH` follows `HIGHSINT64`), `real_t` is
//! `f32` (`REALTYPEWIDTH 32`): every `real_t` expression keeps C's usual
//! arithmetic conversions.
//!
//! Ownership: a `graph_t` owns its arrays (`Vec`); the caller's const
//! `xadj` / `adjncy` that `SetupGraph` stores in `xadj_orig` / `adjncy_orig`
//! are owned copies, never mutated, read through [`graph_t::xadj_eff`] /
//! [`graph_t::adjncy_eff`] (the C `graph->xadj_orig ? graph->xadj_orig :
//! graph->xadj`). Graph chains (`coarser` / `finer`) live in a [`GraphPool`]
//! arena addressed by [`GraphId`]; functions that walk or create graphs take
//! `pool: &mut GraphPool, graph: GraphId`, the others `graph: &mut graph_t`.
//! The `gk_mcore_t` workspace becomes [`wspace`]: every `iwspacemalloc` is a
//! fresh zeroed `Vec` and `WCOREPUSH` / `WCOREPOP` only track depth.
//! `GK_ERREXIT` is `panic!` (the C aborts); debug printing
//! (`HIGHS_ORDERING_PRINT`, compiled out under `NDEBUG`, and only reached with
//! `dbglvl != 0`) is not ported beyond the todo-stubbed printers.
//!
//! Not ported (memory management replaced by Rust ownership): `gk_mcore*`
//! (`GKlib/mcore.c`), `gk_free` (`GKlib/memory.c`), `?realloc`.

pub mod balance;
pub mod coarsen;
pub mod compress;
pub mod fm;
pub mod gklib;
pub mod graph;
pub mod initpart;
pub mod mcutil;
pub mod mmd;
pub mod ometis;
pub mod options;
pub mod pqueue;
pub mod refine;
pub mod sfm;
pub mod wspace;

pub use gklib::{ikv_t, rkv_t};
pub use pqueue::rpq_t;

/// HiGHS: metis.h `idx_t` (`IDXTYPEWIDTH` = 32, or 64 under `HIGHSINT64`)
pub type idx_t = crate::HighsInt;
/// HiGHS: metis.h `real_t` (`REALTYPEWIDTH 32`)
pub type real_t = f32;

/// HiGHS: metis.h `IDX_MAX`
pub const IDX_MAX: idx_t = idx_t::MAX;
/// HiGHS: metis.h `REAL_MAX` (`FLT_MAX`)
pub const REAL_MAX: real_t = f32::MAX;

/// HiGHS: metis.h `METIS_NOPTIONS`
pub const METIS_NOPTIONS: usize = 13;

// HiGHS: metis.h `rstatus_et`
pub const METIS_OK: i32 = 1;
pub const METIS_ERROR_INPUT: i32 = -2;
pub const METIS_ERROR_MEMORY: i32 = -3;
pub const METIS_ERROR: i32 = -4;

// HiGHS: metis.h `moptype_et`
pub const METIS_OP_OMETIS: i32 = 0;

// HiGHS: metis.h `moptions_et`
pub const METIS_OPTION_CTYPE: usize = 0;
pub const METIS_OPTION_IPTYPE: usize = 1;
pub const METIS_OPTION_RTYPE: usize = 2;
pub const METIS_OPTION_DBGLVL: usize = 3;
pub const METIS_OPTION_NITER: usize = 4;
pub const METIS_OPTION_SEED: usize = 5;
pub const METIS_OPTION_COMPRESS: usize = 6;
pub const METIS_OPTION_CCORDER: usize = 7;
pub const METIS_OPTION_PFACTOR: usize = 8;
pub const METIS_OPTION_NSEPS: usize = 9;
pub const METIS_OPTION_UFACTOR: usize = 10;
pub const METIS_OPTION_DROPEDGES: usize = 11;
pub const METIS_OPTION_NO2HOP: usize = 12;

// HiGHS: metis.h `mctype_et`, `miptype_et`, `mrtype_et`
pub const METIS_CTYPE_RM: idx_t = 0;
pub const METIS_CTYPE_SHEM: idx_t = 1;
pub const METIS_IPTYPE_EDGE: idx_t = 0;
pub const METIS_IPTYPE_NODE: idx_t = 1;
pub const METIS_RTYPE_SEP2SIDED: idx_t = 0;
pub const METIS_RTYPE_SEP1SIDED: idx_t = 1;

// HiGHS: metis.h `mdbglvl_et`
pub const METIS_DBG_INFO: idx_t = 1;
pub const METIS_DBG_COARSEN: idx_t = 2;
pub const METIS_DBG_REFINE: idx_t = 4;
pub const METIS_DBG_IPART: idx_t = 8;
pub const METIS_DBG_MOVEINFO: idx_t = 16;
pub const METIS_DBG_SEPINFO: idx_t = 32;

// HiGHS: libmetis/defs.h
pub const HTLENGTH: idx_t = (1 << 13) - 1;
pub const INIT_MAXNAD: idx_t = 200;
pub const UNMATCHED: idx_t = -1;
pub const LARGENIPARTS: idx_t = 7;
pub const SMALLNIPARTS: idx_t = 5;
pub const COARSEN_FRACTION: f64 = 0.85;
pub const COMPRESSION_FRACTION: f64 = 0.85;
pub const MMDSWITCH: idx_t = 120;
pub const OMETIS_DEFAULT_UFACTOR: idx_t = 200;
/// HiGHS: libmetis/coarsen.c line 14
pub const UNMATCHEDFOR2HOP: f64 = 0.10;

/// HiGHS: libmetis/struct.h `nrinfo_t`
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct nrinfo_t {
    pub edegrees: [idx_t; 2],
}

/// HiGHS: libmetis/struct.h `graph_t` (fields in C order; pointers that are
/// `NULL` in C are empty `Vec`s, `xadj_orig` / `adjncy_orig` are `Option`).
#[derive(Clone, Debug, Default)]
pub struct graph_t {
    pub nvtxs: idx_t,
    pub nedges: idx_t,
    pub ncon: idx_t,
    pub xadj: Vec<idx_t>,
    pub vwgt: Vec<idx_t>,
    pub vsize: Vec<idx_t>,
    pub adjncy: Vec<idx_t>,
    pub adjwgt: Vec<idx_t>,
    pub xadj_orig: Option<Vec<idx_t>>,
    pub adjncy_orig: Option<Vec<idx_t>>,
    pub tvwgt: Vec<idx_t>,
    pub invtvwgt: Vec<real_t>,
    pub free_xadj: i32,
    pub free_vwgt: i32,
    pub free_vsize: i32,
    pub free_adjncy: i32,
    pub free_adjwgt: i32,
    pub cmap: Vec<idx_t>,
    pub label: Vec<idx_t>,
    pub mincut: idx_t,
    pub minvol: idx_t,
    pub where_: Vec<idx_t>,
    pub pwgts: Vec<idx_t>,
    pub nbnd: idx_t,
    pub bndptr: Vec<idx_t>,
    pub bndind: Vec<idx_t>,
    pub id: Vec<idx_t>,
    pub ed: Vec<idx_t>,
    pub nrinfo: Vec<nrinfo_t>,
    pub droppedewgt: idx_t,
    pub coarser: Option<GraphId>,
    pub finer: Option<GraphId>,
}

impl graph_t {
    /// `graph->xadj_orig ? graph->xadj_orig : graph->xadj`
    pub fn xadj_eff(&self) -> &[idx_t] {
        match &self.xadj_orig {
            Some(x) => x,
            None => &self.xadj,
        }
    }

    /// `graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy`
    pub fn adjncy_eff(&self) -> &[idx_t] {
        match &self.adjncy_orig {
            Some(a) => a,
            None => &self.adjncy,
        }
    }
}

/// HiGHS: libmetis/struct.h `ctrl_t` (fields in C order; `mcore` is the
/// [`wspace`] depth counter).
#[derive(Clone, Debug, Default)]
pub struct ctrl_t {
    pub optype: i32,
    pub dbglvl: idx_t,
    pub ctype: idx_t,
    pub iptype: idx_t,
    pub rtype: idx_t,
    pub CoarsenTo: idx_t,
    pub no2hop: idx_t,
    pub nseps: idx_t,
    pub ufactor: idx_t,
    pub compress: idx_t,
    pub ccorder: idx_t,
    pub seed: idx_t,
    pub niter: idx_t,
    pub dropedges: idx_t,
    pub maxvwgt: Vec<idx_t>,
    pub ncon: idx_t,
    pub nparts: idx_t,
    pub pfactor: real_t,
    pub ubfactors: Vec<real_t>,
    pub tpwgts: Vec<real_t>,
    pub pijbm: Vec<real_t>,
    pub cfactor: real_t,
    pub mcore: usize,
    pub rng_state: u32,
}

/// Index of a graph in a [`GraphPool`] (the C `graph_t *`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct GraphId(pub usize);

/// Arena holding every live `graph_t` of one `Highs_METIS_NodeND` call.
/// `FreeGraph` empties a slot; slots are never reused, so ids are stable.
#[derive(Clone, Debug, Default)]
pub struct GraphPool {
    pub graphs: Vec<Option<graph_t>>,
}

impl GraphPool {
    pub fn new() -> Self {
        Self::default()
    }

    /// Moves `graph` into the pool (the C `CreateGraph()` result pointer).
    pub fn alloc(&mut self, graph: graph_t) -> GraphId {
        self.graphs.push(Some(graph));
        GraphId(self.graphs.len() - 1)
    }

    pub fn get(&self, id: GraphId) -> &graph_t {
        self.graphs[id.0]
            .as_ref()
            .expect("METIS graph used after FreeGraph")
    }

    pub fn get_mut(&mut self, id: GraphId) -> &mut graph_t {
        self.graphs[id.0]
            .as_mut()
            .expect("METIS graph used after FreeGraph")
    }

    /// Two distinct graphs at once (`ReAdjustMemory(graph, cgraph)`,
    /// `Project2WayNodePartition`).
    pub fn get2_mut(&mut self, a: GraphId, b: GraphId) -> (&mut graph_t, &mut graph_t) {
        assert_ne!(a.0, b.0, "GraphPool::get2_mut needs distinct graphs");
        if a.0 < b.0 {
            let (lo, hi) = self.graphs.split_at_mut(b.0);
            (
                lo[a.0].as_mut().expect("freed"),
                hi[0].as_mut().expect("freed"),
            )
        } else {
            let (lo, hi) = self.graphs.split_at_mut(a.0);
            let (ga, gb) = (
                hi[0].as_mut().expect("freed"),
                lo[b.0].as_mut().expect("freed"),
            );
            (ga, gb)
        }
    }

    /// Drops the graph (`gk_free` of every array plus the struct).
    pub fn release(&mut self, id: GraphId) -> Option<graph_t> {
        self.graphs[id.0].take()
    }

    pub fn is_live(&self, id: GraphId) -> bool {
        self.graphs.get(id.0).is_some_and(|g| g.is_some())
    }
}

/// HiGHS: libmetis/macros.h `GETOPTION`
pub fn GETOPTION(options: Option<&[idx_t]>, idx: usize, defval: idx_t) -> idx_t {
    match options {
        Some(o) if o[idx] != -1 => o[idx],
        _ => defval,
    }
}

/// HiGHS: libmetis/macros.h `I2RUBFACTOR` (a `double` expression)
pub fn I2RUBFACTOR(ufactor: idx_t) -> f64 {
    1.0 + 0.001 * (ufactor as f64)
}

/// HiGHS: libmetis/macros.h::ListInsert (lines 36-40)
/// `lind[n] = i; lptr[i] = n++;` (`BNDInsert` is this macro).
pub fn ListInsert(n: &mut idx_t, lind: &mut [idx_t], lptr: &mut [idx_t], i: idx_t) {
    lind[*n as usize] = i;
    lptr[i as usize] = *n;
    *n += 1;
}

/// HiGHS: libmetis/macros.h::ListDelete (lines 42-47)
/// `lind[lptr[i]] = lind[--n]; lptr[lind[n]] = lptr[i]; lptr[i] = -1;`
/// (`BNDDelete` is this macro).
pub fn ListDelete(n: &mut idx_t, lind: &mut [idx_t], lptr: &mut [idx_t], i: idx_t) {
    *n -= 1;
    lind[lptr[i as usize] as usize] = lind[*n as usize];
    lptr[lind[*n as usize] as usize] = lptr[i as usize];
    lptr[i as usize] = -1;
}

/// HiGHS: libmetis/macros.h `BNDInsert` (lines 53-54)
pub fn BNDInsert(nbnd: &mut idx_t, bndind: &mut [idx_t], bndptr: &mut [idx_t], vtx: idx_t) {
    ListInsert(nbnd, bndind, bndptr, vtx)
}

/// HiGHS: libmetis/macros.h `BNDDelete` (lines 56-57)
pub fn BNDDelete(nbnd: &mut idx_t, bndind: &mut [idx_t], bndptr: &mut [idx_t], vtx: idx_t) {
    ListDelete(nbnd, bndind, bndptr, vtx)
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn list_insert_delete_roundtrip() {
        let mut n: idx_t = 0;
        let mut lind = vec![0; 4];
        let mut lptr = vec![-1; 4];
        ListInsert(&mut n, &mut lind, &mut lptr, 2);
        ListInsert(&mut n, &mut lind, &mut lptr, 0);
        assert_eq!((n, lind[0], lind[1], lptr[2], lptr[0]), (2, 2, 0, 0, 1));
        ListDelete(&mut n, &mut lind, &mut lptr, 2);
        assert_eq!((n, lind[0], lptr[0], lptr[2]), (1, 0, 0, -1));
    }
}
