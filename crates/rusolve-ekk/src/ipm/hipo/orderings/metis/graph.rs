//! HiGHS source: extern/metis/libmetis/graph.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Graph construction and release. `SetupGraph` stores owned copies of the
//! caller's `xadj` / `adjncy` in `xadj_orig` / `adjncy_orig` (the C keeps the
//! const pointers and leaves `xadj` / `adjncy` NULL).

use super::gklib::{imalloc, ismalloc, isum, rmalloc};
use super::{GraphId, GraphPool, ctrl_t, graph_t, idx_t, real_t};

/// HiGHS: graph.c::SetupGraph (lines 17-77)
/// `vwgt` / `vsize` / `adjwgt` are `None` for the C `NULL` (HiPO passes
/// `nullptr` for `vwgt`).
pub fn SetupGraph(
    ctrl: &mut ctrl_t,
    nvtxs: idx_t,
    ncon: idx_t,
    xadj: &[idx_t],
    adjncy: &[idx_t],
    vwgt: Option<&[idx_t]>,
    vsize: Option<&[idx_t]>,
    adjwgt: Option<&[idx_t]>,
) -> graph_t {
    // HiGHS: graph.c lines 20-22 (locals i, j, k, sum, nvwgt, graph)
    let mut graph = CreateGraph();

    // HiGHS: graph.c lines 27-29
    graph.nvtxs = nvtxs;
    graph.nedges = xadj[nvtxs as usize];
    graph.ncon = ncon;

    // HiGHS: graph.c lines 31-33
    graph.xadj_orig = Some(xadj.to_vec());
    graph.xadj = Vec::new();
    graph.free_xadj = 0;

    // HiGHS: graph.c lines 35-37
    graph.adjncy_orig = Some(adjncy.to_vec());
    graph.adjncy = Vec::new();
    graph.free_adjncy = 0;

    // HiGHS: graph.c line 39
    graph.droppedewgt = 0;

    // HiGHS: graph.c lines 42-48 (setup the vertex weights)
    // In Rust, vwgt is Option<&[idx_t]>; if None, we allocate owned.
    // The C code reassigns the local vwgt pointer; we track it via graph.vwgt.
    if let Some(vw) = vwgt {
        graph.vwgt = vw.to_vec();
        graph.free_vwgt = 0;
    } else {
        // C: vwgt = graph->vwgt = ismalloc(ncon*nvtxs, 1);
        graph.vwgt = ismalloc((ncon * nvtxs) as usize, 1);
        graph.free_vwgt = 1;
    }

    // HiGHS: graph.c lines 50-55
    graph.tvwgt = imalloc(ncon as usize);
    graph.invtvwgt = rmalloc(ncon as usize);
    for i in 0..ncon {
        // C: graph->tvwgt[i] = isum(nvtxs, vwgt+i, ncon);
        // vwgt+i means the slice starting at index i, with stride ncon.
        // isum(n, x, incx) sums x[0], x[incx], x[2*incx], ... for n elements.
        // We need to pass a slice starting at i. Since isum takes &[idx_t],
        // we need to be careful. The gklib isum signature is:
        // pub fn isum(n: usize, x: &[idx_t], incx: usize) -> idx_t
        // It accesses x[j * incx] for j in 0..n.
        // So we pass &graph.vwgt[i..] and stride ncon.
        graph.tvwgt[i as usize] = isum(nvtxs as usize, &graph.vwgt[i as usize..], ncon as usize);
        // C: graph->invtvwgt[i] = 1.0/(graph->tvwgt[i] > 0 ? graph->tvwgt[i] : 1);
        // real_t is f32. 1.0 is a double literal in C, but assigned to real_t (f32).
        // The division: 1.0 / (tvwgt[i] > 0 ? tvwgt[i] : 1)
        // tvwgt[i] is idx_t (i32). The ternary yields i32.
        // 1.0 / i32 promotes i32 to double, divides, then truncates to f32.
        let denom: f64 = if graph.tvwgt[i as usize] > 0 {
            graph.tvwgt[i as usize] as f64
        } else {
            1.0
        };
        graph.invtvwgt[i as usize] = (1.0 / denom) as real_t;
    }

    // HiGHS: graph.c lines 59-65 (setup the edge weights)
    if let Some(aw) = adjwgt {
        graph.adjwgt = aw.to_vec();
        graph.free_adjwgt = 0;
    } else {
        // C: adjwgt = graph->adjwgt = ismalloc(graph->nedges, 1);
        graph.adjwgt = ismalloc(graph.nedges as usize, 1);
        graph.free_adjwgt = 1;
    }

    // HiGHS: graph.c line 70
    SetupGraph_tvwgt(&mut graph);

    // HiGHS: graph.c lines 72-73
    if ctrl.optype == crate::ipm::hipo::orderings::metis::METIS_OP_OMETIS {
        SetupGraph_label(&mut graph);
    }

    // HiGHS: graph.c line 76
    graph
}

/// HiGHS: graph.c::SetupGraph_tvwgt (lines 83-96)
pub fn SetupGraph_tvwgt(graph: &mut graph_t) {
    // HiGHS: graph.c lines 87-90
    if graph.tvwgt.is_empty() {
        graph.tvwgt = imalloc(graph.ncon as usize);
    }
    if graph.invtvwgt.is_empty() {
        graph.invtvwgt = rmalloc(graph.ncon as usize);
    }

    // HiGHS: graph.c lines 92-95
    for i in 0..graph.ncon {
        // C: graph->tvwgt[i] = isum(graph->nvtxs, graph->vwgt+i, graph->ncon);
        graph.tvwgt[i as usize] = isum(
            graph.nvtxs as usize,
            &graph.vwgt[i as usize..],
            graph.ncon as usize,
        );
        // C: graph->invtvwgt[i] = 1.0/(graph->tvwgt[i] > 0 ? graph->tvwgt[i] : 1);
        let denom: f64 = if graph.tvwgt[i as usize] > 0 {
            graph.tvwgt[i as usize] as f64
        } else {
            1.0
        };
        graph.invtvwgt[i as usize] = (1.0 / denom) as real_t;
    }
}

/// HiGHS: graph.c::SetupGraph_label (lines 102-111)
pub fn SetupGraph_label(graph: &mut graph_t) {
    // HiGHS: graph.c lines 106-107
    if graph.label.is_empty() {
        graph.label = imalloc(graph.nvtxs as usize);
    }

    // HiGHS: graph.c lines 109-110
    for i in 0..graph.nvtxs {
        graph.label[i as usize] = i;
    }
}

/// HiGHS: graph.c::SetupSplitGraph (lines 117-140)
pub fn SetupSplitGraph(graph: &graph_t, snvtxs: idx_t, snedges: idx_t) -> graph_t {
    // HiGHS: graph.c line 121
    let mut sgraph = CreateGraph();

    // HiGHS: graph.c lines 123-125
    sgraph.nvtxs = snvtxs;
    sgraph.nedges = snedges;
    sgraph.ncon = graph.ncon;

    // HiGHS: graph.c lines 128-134
    sgraph.xadj = imalloc((snvtxs + 1) as usize);
    sgraph.vwgt = imalloc((sgraph.ncon * snvtxs) as usize);
    sgraph.adjncy = imalloc(snedges as usize);
    sgraph.adjwgt = imalloc(snedges as usize);
    sgraph.label = imalloc(snvtxs as usize);
    sgraph.tvwgt = imalloc(sgraph.ncon as usize);
    sgraph.invtvwgt = rmalloc(sgraph.ncon as usize);

    // HiGHS: graph.c lines 136-137
    if !graph.vsize.is_empty() {
        sgraph.vsize = imalloc(snvtxs as usize);
    }

    // HiGHS: graph.c line 139
    sgraph
}

/// HiGHS: graph.c::CreateGraph (lines 146-155)
pub fn CreateGraph() -> graph_t {
    // HiGHS: graph.c lines 150-152
    // C: malloc(sizeof(graph_t)) then InitGraph(graph)
    // In Rust, we create a default graph_t and then call InitGraph to set
    // the sentinel values. Default gives zeros/empty vecs; InitGraph overwrites
    // the scalar fields with -1 and the free_* flags with 1.
    let mut graph = graph_t::default();
    InitGraph(&mut graph);
    graph
}

/// HiGHS: graph.c::InitGraph (lines 161-208)
pub fn InitGraph(graph: &mut graph_t) {
    // HiGHS: graph.c line 163 (memset to 0) — Rust Default already zeroes.
    // We still set the sentinel values explicitly.

    // HiGHS: graph.c lines 166-171 (graph size constants)
    graph.nvtxs = -1;
    graph.nedges = -1;
    graph.ncon = -1;
    graph.mincut = -1;
    graph.minvol = -1;
    graph.nbnd = -1;

    // HiGHS: graph.c lines 174-182 (memory for the graph structure)
    // C sets pointers to NULL; in Rust, empty Vec is the equivalent.
    graph.xadj = Vec::new();
    graph.vwgt = Vec::new();
    graph.vsize = Vec::new();
    graph.adjncy = Vec::new();
    graph.adjwgt = Vec::new();
    graph.label = Vec::new();
    graph.cmap = Vec::new();
    graph.tvwgt = Vec::new();
    graph.invtvwgt = Vec::new();

    // HiGHS: graph.c lines 184-185
    graph.xadj_orig = None;
    graph.adjncy_orig = None;

    // HiGHS: graph.c lines 188-192 (free flags default to 1)
    graph.free_xadj = 1;
    graph.free_vwgt = 1;
    graph.free_vsize = 1;
    graph.free_adjncy = 1;
    graph.free_adjwgt = 1;

    // HiGHS: graph.c lines 196-202 (partition/refinement structure)
    graph.where_ = Vec::new();
    graph.pwgts = Vec::new();
    graph.id = Vec::new();
    graph.ed = Vec::new();
    graph.bndptr = Vec::new();
    graph.bndind = Vec::new();
    graph.nrinfo = Vec::new();

    // HiGHS: graph.c lines 205-206 (linked-list structure)
    graph.coarser = None;
    graph.finer = None;
}

/// HiGHS: graph.c::FreeSData (lines 214-227)
pub fn FreeSData(graph: &mut graph_t) {
    // HiGHS: graph.c lines 217-218
    if graph.free_xadj != 0 {
        graph.xadj = Vec::new();
    }
    // HiGHS: graph.c lines 219-220
    if graph.free_vwgt != 0 {
        graph.vwgt = Vec::new();
    }
    // HiGHS: graph.c lines 221-222
    if graph.free_vsize != 0 {
        graph.vsize = Vec::new();
    }
    // HiGHS: graph.c lines 223-224
    if graph.free_adjncy != 0 {
        graph.adjncy = Vec::new();
    }
    // HiGHS: graph.c lines 225-226
    if graph.free_adjwgt != 0 {
        graph.adjwgt = Vec::new();
    }
}

/// HiGHS: graph.c::FreeRData (lines 233-243)
pub fn FreeRData(graph: &mut graph_t) {
    // HiGHS: graph.c lines 236-242
    graph.where_ = Vec::new();
    graph.pwgts = Vec::new();
    graph.id = Vec::new();
    graph.ed = Vec::new();
    graph.bndptr = Vec::new();
    graph.bndind = Vec::new();
    graph.nrinfo = Vec::new();
}

/// HiGHS: graph.c::FreeGraph (lines 249-268)
/// Releases the pool slot and sets `*r_graph = NULL`.
pub fn FreeGraph(pool: &mut GraphPool, r_graph: &mut Option<GraphId>) {
    // HiGHS: graph.c line 253
    let gid = r_graph.take().expect("FreeGraph called on NULL graph");

    // HiGHS: graph.c lines 256-259
    // We need to get the graph, free its data, then release the slot.
    // FreeSData and FreeRData operate on the graph in place.
    {
        let graph = pool.get_mut(gid);
        FreeSData(graph);
        FreeRData(graph);
    }

    // HiGHS: graph.c lines 261-264
    {
        let graph = pool.get_mut(gid);
        graph.tvwgt = Vec::new();
        graph.invtvwgt = Vec::new();
        graph.label = Vec::new();
        graph.cmap = Vec::new();
    }

    // HiGHS: graph.c line 265 (gk_free(&graph))
    pool.release(gid);

    // HiGHS: graph.c line 267 (*r_graph = NULL)
    // Already done by take() above.
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn setup_graph_unit_weights() {
        let mut ctrl = ctrl_t {
            optype: crate::ipm::hipo::orderings::metis::METIS_OP_OMETIS,
            ..Default::default()
        };
        // path 0-1-2
        let xadj = [0, 1, 3, 4];
        let adjncy = [1, 0, 2, 1];
        let g = SetupGraph(&mut ctrl, 3, 1, &xadj, &adjncy, None, None, None);
        assert_eq!((g.nvtxs, g.nedges, g.tvwgt[0]), (3, 4, 3));
        assert_eq!(g.label, vec![0, 1, 2]);
        assert_eq!(g.xadj_eff(), &xadj);
    }
}
