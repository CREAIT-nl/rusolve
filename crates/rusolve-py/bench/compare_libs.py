"""Time how long it takes to *build* the same LP in Python through different interfaces.

Model (facility-location LP relaxation), F facilities x C customers:
  vars  y[f] in [0,1] cost fixed[f];  x[f,c] in [0,1] cost cost[f,c]
  rows  sum_f x[f,c] == 1                (C rows, F nnz each)
        x[f,c] - y[f] <= 0               (F*C rows, 2 nnz each)
Column order: y first (F), then x row-major (f*C + c).
"""
import sys, time, gc
import numpy as np

F, C = (int(a) for a in sys.argv[1:3]) if len(sys.argv) >= 3 else (100, 1000)
ONLY = set(sys.argv[3].split(",")) if len(sys.argv) >= 4 else None
rng = np.random.default_rng(0)
fixed = rng.uniform(10, 100, F)
cost = rng.uniform(1, 20, (F, C))
NCOL, NROW, NNZ = F + F * C, C + F * C, F * C + 2 * F * C
INF = 1e30


def csr_numpy():
    """Vectorised CSR assembly: what a user (or a fast layer) has to compute anyway."""
    xidx = F + np.arange(F * C).reshape(F, C)          # column index of x[f,c]
    # demand rows: row c has x[0..F-1, c]
    d_idx = xidx.T.ravel()
    d_val = np.ones(F * C)
    d_start = np.arange(0, F * C + 1, F)
    # link rows: row (f,c) has x[f,c] (+1), y[f] (-1)
    l_idx = np.empty(2 * F * C, dtype=np.int32)
    l_idx[0::2] = xidx.ravel()
    l_idx[1::2] = np.repeat(np.arange(F), C)
    l_val = np.empty(2 * F * C)
    l_val[0::2], l_val[1::2] = 1.0, -1.0
    start = np.concatenate([d_start, F * C + np.arange(2, 2 * F * C + 1, 2)]).astype(np.int32)
    index = np.concatenate([d_idx, l_idx]).astype(np.int32)
    value = np.concatenate([d_val, l_val])
    rlo = np.concatenate([np.ones(C), np.full(F * C, -INF)])
    rup = np.concatenate([np.ones(C), np.zeros(F * C)])
    ccost = np.concatenate([fixed, cost.ravel()])
    clo, cup = np.zeros(NCOL), np.ones(NCOL)
    return ccost, clo, cup, rlo, rup, start, index, value


# ---------------------------------------------------------------- highspy
def hs_new():
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    return h


def hs_expr_loop():
    """Idiomatic highspy modelling API, one addConstr per row."""
    h = hs_new()
    y = h.addVariables(F, lb=0, ub=1, obj=fixed.tolist())
    x = h.addVariables(F, C, lb=0, ub=1, obj=cost.ravel().tolist())
    for c in range(C):
        h.addConstr(h.qsum(x[f, c] for f in range(F)) == 1)
    for f in range(F):
        for c in range(C):
            h.addConstr(x[f, c] <= y[f])
    return h


def hs_expr_vectorised():
    """highspy HighspyArray (numpy object arrays of expressions)."""
    h = hs_new()
    y = h.addVariables(F, lb=0, ub=1, obj=fixed.tolist())
    x = h.addVariables(F, C, lb=0, ub=1, obj=cost.ravel().tolist())
    h.addConstrs(x.sum(axis=0) == 1)
    h.addConstrs((x <= y[:, None]).ravel())
    return h


def hs_addrow_loop():
    """Arrays precomputed; one pybind11 addRow call per row (boundary cost only)."""
    import highspy
    ccost, clo, cup, rlo, rup, start, index, value = csr_numpy()
    h = hs_new()
    h.addCols(NCOL, ccost, clo, cup, 0, np.array([], np.int32), np.array([], np.int32), np.array([]))
    for r in range(NROW):
        s, e = start[r], start[r + 1]
        h.addRow(rlo[r], rup[r], e - s, index[s:e], value[s:e])
    return h


def hs_addrows_batch():
    """numpy CSR, one addCols + one addRows call."""
    ccost, clo, cup, rlo, rup, start, index, value = csr_numpy()
    h = hs_new()
    h.addCols(NCOL, ccost, clo, cup, 0, np.array([], np.int32), np.array([], np.int32), np.array([]))
    h.addRows(NROW, rlo, rup, len(index), start[:-1], index, value)
    return h


def hs_passmodel_flat():
    import highspy
    ccost, clo, cup, rlo, rup, start, index, value = csr_numpy()
    h = hs_new()
    h.passModel(NCOL, NROW, len(index), int(highspy.MatrixFormat.kRowwise), int(highspy.ObjSense.kMinimize),
                0.0, ccost, clo, cup, rlo, rup, start, index, value, np.zeros(NCOL, np.int32))
    return h


def hs_highslp():
    import highspy
    ccost, clo, cup, rlo, rup, start, index, value = csr_numpy()
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = NCOL, NROW
    lp.col_cost_, lp.col_lower_, lp.col_upper_ = ccost, clo, cup
    lp.row_lower_, lp.row_upper_ = rlo, rup
    lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_, lp.a_matrix_.index_, lp.a_matrix_.value_ = start, index, value
    h = hs_new()
    h.passModel(lp)
    return h


# ---------------------------------------------------------------- rusolve-py (current crate)
def rs_new():
    import rusolve
    h = rusolve.Highs()
    h.setOptionValue("output_flag", False)
    return h


def rs_flat_numpy():
    import rusolve
    ccost, clo, cup, rlo, rup, start, index, value = csr_numpy()
    h = rs_new()
    h.passModel(NCOL, NROW, len(index), rusolve.MatrixFormat.kRowwise, rusolve.ObjSense.kMinimize,
                0.0, ccost, clo, cup, rlo, rup, start, index, value)
    return h


def rs_flat_list():
    import rusolve
    a = [v.tolist() for v in csr_numpy()]
    h = rs_new()
    h.passModel(NCOL, NROW, len(a[6]), rusolve.MatrixFormat.kRowwise, rusolve.ObjSense.kMinimize, 0.0, *a)
    return h


def rs_addrow_loop():
    ccost, clo, cup, rlo, rup, start, index, value = csr_numpy()
    h = rs_new()
    for j in range(NCOL):
        h.addCol(float(ccost[j]), 0.0, 1.0, 0, [], [])
    idx, val = index.tolist(), value.tolist()
    for r in range(NROW):
        s, e = int(start[r]), int(start[r + 1])
        h.addRow(float(rlo[r]), float(rup[r]), e - s, idx[s:e], val[s:e])
    return h


# ---------------------------------------------------------------- other modelling layers
def poi_loop():
    """pyoptinterface: C++ expression objects, one call per row."""
    import pyoptinterface as poi
    from pyoptinterface import highs
    highs.autoload_library()
    m = highs.Model()
    m.set_model_attribute(poi.ModelAttribute.Silent, True)
    y = [m.add_variable(lb=0, ub=1) for f in range(F)]
    x = [[m.add_variable(lb=0, ub=1) for c in range(C)] for f in range(F)]
    obj = poi.ExprBuilder()
    for f in range(F):
        obj.add_affine_term(y[f], fixed[f])
        for c in range(C):
            obj.add_affine_term(x[f][c], cost[f, c])
    m.set_objective(obj, poi.ObjectiveSense.Minimize)
    for c in range(C):
        m.add_linear_constraint(poi.quicksum(x[f][c] for f in range(F)), poi.Eq, 1.0)
    for f in range(F):
        yf = y[f]
        for c in range(C):
            m.add_linear_constraint(x[f][c] - yf, poi.Leq, 0.0)
    return m


def linopy_build():
    """linopy: xarray-vectorised; build + conversion to a highspy object (no solve)."""
    import linopy, xarray as xr, pandas as pd
    m = linopy.Model()
    fi, ci = xr.DataArray(range(F), dims="f"), xr.DataArray(range(C), dims="c")
    y = m.add_variables(0, 1, coords=[pd.Index(range(F), name='f')], name="y")
    x = m.add_variables(0, 1, coords=[pd.Index(range(F), name='f'), pd.Index(range(C), name='c')], name="x")
    m.add_constraints(x.sum("f") == 1, name="demand")
    m.add_constraints(x - y <= 0, name="link")
    m.add_objective((xr.DataArray(fixed, coords=[fi]) * y).sum()
                    + (xr.DataArray(cost, coords=[fi, ci]) * x).sum())
    t = time.perf_counter()
    h = m.to_highspy()
    linopy_build.convert = time.perf_counter() - t
    return h


def pulp_build():
    import pulp
    p = pulp.LpProblem("fl", pulp.LpMinimize)
    y = [p.add_variable(f"y{f}", lowBound=0, upBound=1) for f in range(F)]
    x = [[p.add_variable(f"x{f}_{c}", lowBound=0, upBound=1) for c in range(C)] for f in range(F)]
    p += pulp.lpSum(fixed[f] * y[f] for f in range(F)) + pulp.lpSum(
        cost[f, c] * x[f][c] for f in range(F) for c in range(C))
    for c in range(C):
        p += pulp.lpSum(x[f][c] for f in range(F)) == 1
    for f in range(F):
        for c in range(C):
            p += x[f][c] <= y[f]
    import highspy
    t = time.perf_counter()
    p.solverModel = highspy.Highs()
    pulp.HiGHS(msg=False).buildSolverModel(p)
    pulp_build.convert = time.perf_counter() - t
    return p


def dims(m):
    for get in (lambda: (m.getNumCol(), m.getNumRow(), m.getLp().a_matrix_.numNz()),
                lambda: (m.getNumCol(), m.getNumRow(), len(m.getLp().a_matrix_.index_)),
                lambda: (m.number_of_variables(), m.number_of_constraints(poi_eq()), None),
                lambda: (len(m.variables()), m.numConstraints(), None)):
        try:
            return get()
        except Exception:
            pass
    return None


def poi_eq():
    import pyoptinterface as poi
    return poi.ConstraintType.Linear


def poi_le():
    return poi_eq()


METHODS = [csr_numpy, hs_passmodel_flat, hs_highslp, hs_addrows_batch, hs_addrow_loop,
           hs_expr_vectorised, hs_expr_loop, rs_flat_numpy, rs_flat_list, rs_addrow_loop,
           poi_loop, linopy_build, pulp_build]

if __name__ == "__main__":
    print(f"F={F} C={C}: {NCOL} cols, {NROW} rows, {NNZ} nnz")
    for fn in METHODS:
        if ONLY and fn.__name__ not in ONLY:
            continue
        gc.collect()
        t = time.perf_counter()
        try:
            m = fn()
        except Exception as e:
            print(f"{fn.__name__:22s} FAILED {type(e).__name__}: {e}")
            continue
        dt = time.perf_counter() - t
        extra = f"  (of which conversion to highspy {fn.convert:.3f}s)" if hasattr(fn, "convert") else ""
        d = dims(m) if fn is not csr_numpy else ""
        print(f"{fn.__name__:22s} {dt:9.3f} s  {1e9 * dt / NNZ:8.1f} ns/nnz  {d}{extra}", flush=True)
        del m
