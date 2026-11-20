import time, numpy as np, scipy.sparse as sp, sys
sys.argv = [sys.argv[0], "100", "2000"]
import compare_libs as b, highspy
F, C, NNZ = b.F, b.C, b.NNZ
def t(label, fn, n=NNZ):
    s = time.perf_counter(); fn(); d = time.perf_counter() - s
    print(f"{label:58s} {d:7.3f} s {1e9*d/n:7.1f} ns/nnz")
xs = list(range(F*C)); ys = list(range(F))
def bare_loop():
    for f in range(F):
        yf = ys[f]
        for c in range(C):
            v = xs[f*C+c]
t("floor: bare Python double loop, no objects", bare_loop, F*C*2)
def tuple_terms():
    rows = []
    for f in range(F):
        yf = ys[f]
        for c in range(C):
            rows.append(((xs[f*C+c], 1.0), (yf, -1.0)))
t("floor: loop + build tuple-of-(idx,coef) per row", tuple_terms, F*C*2)
cc, lo, up, rlo, rup, start, index, value = b.csr_numpy()
def addrows_i64():
    h = b.hs_new(); h.addCols(b.NCOL, cc, lo, up, 0, np.array([],np.int32), np.array([],np.int32), np.array([]))
    h.addRows(b.NROW, rlo, rup, len(index), start[:-1].astype(np.int64), index.astype(np.int64), value)
t("highspy addRows with int64 index arrays (incl. astype)", addrows_i64)
A = sp.csr_matrix((value, index, start), shape=(b.NROW, b.NCOL))
def scipy_build():
    # what a matrix-API user writes: blocks via scipy.sparse
    Ix = sp.identity(F*C, format="csr")
    Y = sp.kron(sp.identity(F), np.ones((C,1)), format="csr")
    D = sp.kron(np.ones((1,F)), sp.identity(C), format="csr")
    M = sp.vstack([sp.hstack([sp.csr_matrix((C,F)), D]), sp.hstack([-Y, Ix])], format="csr")
    return M
t("user side: assemble same A with scipy.sparse kron/vstack", scipy_build)
def solver_side_copy():
    a = np.array(value, copy=True); i = np.array(index, copy=True); s_ = np.array(start, copy=True)
t("floor: memcpy of CSR arrays (what zero-copy ingest costs)", solver_side_copy)
