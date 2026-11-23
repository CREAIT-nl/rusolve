"""Acceptance benchmark for the `rusolve.Model` API: targets P1-P6 of docs/design/python-api.md.

    python bench_api.py                 # all, full size (P5 builds a 3M-nnz model)
    python bench_api.py --quick         # short run for CI-like use (small sizes, few reps)
    python bench_api.py p1 p4           # only some of them; names: p1 p2 p3 p4 p5 p6
    python bench_api.py --no-compare    # skip highspy / pyoptinterface rows

Facility-location model (same as compare_libs.py), F facilities x C customers:
  y[f] in [0,1] cost fixed[f]; x[f,c] in [0,1] cost cost[f,c]
  sum_f x[f,c] == 1 (C rows, F nnz each);  x[f,c] - y[f] <= 0 (F*C rows, 2 nnz each)
  NNZ = 3 F C.  Column order: y first, then x row-major.

Every row is "PASS"/"FAIL" against its target; the exit status is 1 if any FAIL.
Quick mode runs the same checks at small sizes, where fixed costs dominate: its P3/P4/P5 values
are indicative only.
"""
import argparse
import gc
import statistics
import sys
import time

import numpy as np
import rusolve as rs

INF = float("inf")
NS = 1e-9
RESULTS = []  # (id, description, value, unit, target, ok)


# ------------------------------------------------------------------------------- reporting
def report(pid, what, value, unit, target=None, limit="le"):
    ok = None
    tgt = ""
    if target is not None:
        ok = value <= target if limit == "le" else value >= target
        tgt = f"{'<=' if limit == 'le' else '>='} {target:g} {unit}"
    mark = "" if ok is None else ("PASS" if ok else "FAIL")
    num = f"{value:12.2f}" if unit in ("x", "ms") else f"{value:12.1f}"
    print(f"{pid:3s} {what:58s} {num} {unit:8s} {tgt:20s} {mark}", flush=True)
    if ok is not None:
        RESULTS.append((pid, what, value, unit, target, ok))


def header(title):
    print(f"\n== {title}", flush=True)


def median_time_ns(fn, reps, warmup=20):
    for _ in range(warmup):
        fn()
    ts = []
    gc.collect()
    gc.disable()
    try:
        for _ in range(reps):
            t = time.perf_counter_ns()
            fn()
            ts.append(time.perf_counter_ns() - t)
    finally:
        gc.enable()
    return statistics.median(ts)


def timed(fn, reps):
    """Median wall time in seconds of fn() over reps runs, each run followed by a collect."""
    ts, out = [], None
    for _ in range(reps):
        gc.collect()
        t = time.perf_counter()
        out = fn()
        ts.append(time.perf_counter() - t)
        del out
    return statistics.median(ts)


def force_build(m):
    """Make the Model hand its HighsLp to the solver, without searching.

    Uses a dedicated method if the binding exposes one; otherwise solves with a zero time limit.
    Returns the name of what was used."""
    for name in ("_sync", "build", "_build", "finalize", "_finalize", "to_lp"):
        f = getattr(m, name, None)
        if callable(f):
            f()
            return name
    try:
        m.solve(time_limit=0.0)
    except (ValueError, RuntimeError):
        m.solve(time_limit=1e-9)
    return "solve(time_limit=0)"


# ------------------------------------------------------------------------------------ P1
def p1(args):
    header("P1  tiny model round trip: 1 var, 1 row, build + solve + read objective (median)")
    reps = 300 if args.quick else 3000

    def rs_model():
        m = rs.Model()
        x = m.add_var(cost=-1.0)
        m.add_constr({x: 1.0}, "<=", 3.0)
        return m.solve().objective

    assert abs(rs_model() + 3.0) < 1e-9
    t_rs = median_time_ns(rs_model, reps) / 1e3
    report("P1", "rusolve.Model", t_rs, "us", 20.0)

    import rusolve.scipy as rss
    cvec = [-1.0]
    cons_args = ([[1.0]], -INF, 3.0)

    class LC:  # a duck-typed LinearConstraint, so the shim is timed without scipy objects
        def __init__(self, A, lb, ub):
            self.A, self.lb, self.ub = A, lb, ub

    def shim():
        return rss.milp(cvec, constraints=LC(*cons_args)).fun

    assert abs(shim() + 3.0) < 1e-9
    t_shim = median_time_ns(shim, reps) / 1e3
    report("P1", "rusolve.scipy.milp (duck-typed constraint)", t_shim, "us")
    report("P1", "  shim overhead over rusolve.Model", t_shim - t_rs, "us")

    if args.compare:
        try:
            from scipy.optimize import LinearConstraint, milp

            def sp_milp():
                return milp(cvec, constraints=LinearConstraint([[1.0]], -INF, 3.0)).fun

            report("P1", "scipy.optimize.milp (for scale)",
                   median_time_ns(sp_milp, max(reps // 3, 100)) / 1e3, "us")

            def shim_sp():
                return rss.milp(cvec, constraints=LinearConstraint([[1.0]], -INF, 3.0)).fun

            report("P1", "rusolve.scipy.milp (scipy LinearConstraint)",
                   median_time_ns(shim_sp, reps) / 1e3, "us")
        except ImportError:
            print("P1  scipy missing, row skipped")
        try:
            import highspy

            def hs():
                h = highspy.Highs()
                h.setOptionValue("output_flag", False)
                h.addVar(0.0, INF)
                h.changeColCost(0, -1.0)
                h.addRow(-INF, 3.0, 1, [0], [1.0])
                h.run()
                return h.getInfo().objective_function_value

            assert abs(hs() + 3.0) < 1e-9
            report("P1", "highspy (for scale)", median_time_ns(hs, max(reps // 3, 100)) / 1e3, "us")
        except ImportError:
            print("P1  highspy missing, row skipped")


# ------------------------------------------------------------------------------------ P2
def p2(args):
    header("P2  add_constrs with dict rows: library time per term")
    n_rows, k = (2000, 20) if args.quick else (20000, 50)
    reps = 3 if args.quick else 5
    ncol = 1000
    rng = np.random.default_rng(0)
    cols = rng.integers(0, ncol, (n_rows, k))
    vals = rng.uniform(1, 2, (n_rows, k))
    nterms = n_rows * k
    gc.collect()
    t0 = time.perf_counter()
    rows = [dict(zip(cols[i].tolist(), vals[i].tolist())) for i in range(n_rows)]
    t_dicts = time.perf_counter() - t0
    nterms = sum(len(r) for r in rows)  # duplicate columns collapse in the dict
    rhs = [10.0] * n_rows

    def add():
        m = rs.Model()
        m.add_vars(ncol)
        t = time.perf_counter()
        m.add_constrs(rows, "<=", rhs)
        return time.perf_counter() - t

    ts = sorted(add() for _ in range(reps))
    t = ts[len(ts) // 2]
    report("P2", f"add_constrs, dict rows ({n_rows} rows x ~{k} terms)", t / nterms / NS, "ns/term",
           30.0)
    report("P2", "  (user cost of building the dicts, not charged)", t_dicts / nterms / NS,
           "ns/term")
    # the same rows as (indices, values) pairs of numpy arrays: no per-term Python objects at all
    pairs = [(cols[i].astype(np.int32), vals[i]) for i in range(n_rows)]

    def add_pairs():
        m = rs.Model()
        m.add_vars(ncol)
        t = time.perf_counter()
        m.add_constrs(pairs, "<=", rhs)
        return time.perf_counter() - t

    ts = sorted(add_pairs() for _ in range(reps))
    report("P2", "  add_constrs, (int32, float64) array pairs", ts[len(ts) // 2] / (n_rows * k) / NS,
           "ns/term")


# ------------------------------------------------------------------------- FL data + builders
class FL:
    def __init__(self, F, C):
        self.F, self.C = F, C
        rng = np.random.default_rng(0)
        self.fixed = rng.uniform(10, 100, F)
        self.cost = rng.uniform(1, 20, (F, C))
        self.ncol, self.nrow, self.nnz = F + F * C, C + F * C, 3 * F * C

    def csr(self):
        """Vectorised CSR (demand rows first, then link rows), as in compare_libs.csr_numpy."""
        F, C = self.F, self.C
        xidx = F + np.arange(F * C).reshape(F, C)
        d_idx = xidx.T.ravel()
        l_idx = np.empty(2 * F * C, dtype=np.int32)
        l_idx[0::2] = xidx.ravel()
        l_idx[1::2] = np.repeat(np.arange(F), C)
        l_val = np.empty(2 * F * C)
        l_val[0::2], l_val[1::2] = 1.0, -1.0
        d_start = np.arange(0, F * C + 1, F)
        l_start = np.arange(0, 2 * F * C + 1, 2)
        return d_start, d_idx.astype(np.int32), np.ones(F * C), l_start, l_idx, l_val


def rs_expr(fl):
    """Expression API: Var arithmetic, quicksum, comparisons."""
    F, C = fl.F, fl.C
    m = rs.Model()
    y = m.add_vars(F, ub=1, cost=fl.fixed)
    x = m.add_vars(F * C, ub=1, cost=fl.cost.ravel())
    for c in range(C):
        m.add_constr(rs.quicksum(x[f * C + c] for f in range(F)) == 1)
    for f in range(F):
        yf = y[f]
        base = f * C
        m.add_constrs(x[base + c] <= yf for c in range(C))
    return m


def rs_dict(fl, rows=None):
    F, C = fl.F, fl.C
    m = rs.Model()
    m.add_vars(F, ub=1, cost=fl.fixed)
    m.add_vars(F * C, ub=1, cost=fl.cost.ravel())
    if rows is None:
        rows = dict_rows(fl)
    dem, link = rows
    m.add_constrs(dem, "==", 1.0)
    m.add_constrs(link, "<=", 0.0)
    return m


def dict_rows(fl):
    F, C = fl.F, fl.C
    dem = [{F + f * C + c: 1.0 for f in range(F)} for c in range(C)]
    link = [{F + f * C + c: 1.0, f: -1.0} for f in range(F) for c in range(C)]
    return dem, link


def rs_csr(fl, parts=None):
    import scipy.sparse as sp
    F, C = fl.F, fl.C
    t0 = time.perf_counter()
    d_start, d_idx, d_val, l_start, l_idx, l_val = parts or fl.csr()
    D = sp.csr_matrix((d_val, d_idx, d_start), shape=(C, fl.ncol))
    L = sp.csr_matrix((l_val, l_idx, l_start), shape=(F * C, fl.ncol))
    rs_csr.assemble = time.perf_counter() - t0
    m = rs.Model()
    m.add_vars(F, ub=1, cost=fl.fixed)
    m.add_vars(F * C, ub=1, cost=fl.cost.ravel())
    m.add_matrix(D, "==", np.ones(C))
    m.add_matrix(L, "<=", np.zeros(F * C))
    return m


def hs_expr(fl):
    import highspy
    F, C = fl.F, fl.C
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    y = h.addVariables(F, lb=0, ub=1, obj=fl.fixed.tolist())
    x = h.addVariables(F, C, lb=0, ub=1, obj=fl.cost.ravel().tolist())
    for c in range(C):
        h.addConstr(h.qsum(x[f, c] for f in range(F)) == 1)
    for f in range(F):
        for c in range(C):
            h.addConstr(x[f, c] <= y[f])
    return h


def hs_csr(fl):
    import highspy
    d_start, d_idx, d_val, l_start, l_idx, l_val = fl.csr()
    start = np.concatenate([d_start, fl.F * fl.C + l_start[1:]]).astype(np.int32)
    index = np.concatenate([d_idx, l_idx]).astype(np.int32)
    value = np.concatenate([d_val, l_val])
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    ccost = np.concatenate([fl.fixed, fl.cost.ravel()])
    z, o = np.zeros(fl.ncol), np.ones(fl.ncol)
    h.addCols(fl.ncol, ccost, z, o, 0, np.array([], np.int32), np.array([], np.int32),
              np.array([]))
    rlo = np.concatenate([np.ones(fl.C), np.full(fl.F * fl.C, -INF)])
    rup = np.concatenate([np.ones(fl.C), np.zeros(fl.F * fl.C)])
    h.addRows(fl.nrow, rlo, rup, len(index), start[:-1], index, value)
    return h


def poi_loop(fl):
    import pyoptinterface as poi
    from pyoptinterface import highs
    F, C = fl.F, fl.C
    highs.autoload_library()
    m = highs.Model()
    m.set_model_attribute(poi.ModelAttribute.Silent, True)
    y = [m.add_variable(lb=0, ub=1) for _ in range(F)]
    x = [[m.add_variable(lb=0, ub=1) for _ in range(C)] for _ in range(F)]
    obj = poi.ExprBuilder()
    for f in range(F):
        obj.add_affine_term(y[f], fl.fixed[f])
        for c in range(C):
            obj.add_affine_term(x[f][c], fl.cost[f, c])
    m.set_objective(obj, poi.ObjectiveSense.Minimize)
    for c in range(C):
        m.add_linear_constraint(poi.quicksum(x[f][c] for f in range(F)), poi.Eq, 1.0)
    for f in range(F):
        for c in range(C):
            m.add_linear_constraint(x[f][c] - y[f], poi.Leq, 0.0)
    return m


def build_rusolve(fn, fl, reps):
    """ns/nnz of fn(fl) including the hand-over of the model to the solver (a built HighsLp)."""
    how = []

    def run():
        m = fn(fl)
        how.append(force_build(m))
        assert m.num_vars == fl.ncol and m.num_constrs == fl.nrow, (m.num_vars, m.num_constrs)
        return m

    t = timed(run, reps)
    return t / fl.nnz / NS, how[0]


def sizes(args):
    if args.quick:
        return (30, 700), (30, 3500)    # ~63k nnz, ~315k nnz
    return (100, 2000), (100, 10000)    # 600k nnz, 3M nnz


# ------------------------------------------------------------------------------------ P3
def p3(args):
    header("P3  add_matrix(scipy csr) into a built model (call to built HighsLp)")
    import scipy.sparse as sp
    F, C = sizes(args)[0]
    fl = FL(F, C)
    parts = fl.csr()
    d_start, d_idx, d_val, l_start, l_idx, l_val = parts
    L = sp.csr_matrix((l_val, l_idx, l_start), shape=(F * C, fl.ncol))
    D = sp.csr_matrix((d_val, d_idx, d_start), shape=(C, fl.ncol))
    nnz = L.nnz + D.nnz
    has_sync = callable(getattr(rs.Model, "_sync", None))
    reps = 3 if args.quick else 5
    t_add, t_build, how = [], [], "none"
    for _ in range(reps):
        gc.collect()
        m = rs.Model()
        m.add_vars(fl.ncol, ub=1)
        t = time.perf_counter()
        m.add_matrix(D, "==", np.ones(C))
        m.add_matrix(L, "<=", np.zeros(F * C))
        t_add.append(time.perf_counter() - t)
        t = time.perf_counter()
        how = force_build(m)
        t_build.append(time.perf_counter() - t)
    a, b = statistics.median(t_add), statistics.median(t_build)
    if has_sync:
        report("P3", f"add_matrix(csr, {nnz} nnz) + build via {how}", (a + b) / nnz / NS, "ns/nnz",
               5.0)
    else:
        # no build-only hook: the hand-over to HiGHS can only be forced by solve(), which also
        # pays HiGHS's own model setup. Target is checked on the add_matrix call alone.
        report("P3", f"add_matrix(csr, {nnz} nnz) call only (no _sync available)", a / nnz / NS,
               "ns/nnz", 5.0)
        report("P3", f"  hand-over forced via {how} (includes HiGHS setup)", b / nnz / NS,
               "ns/nnz")


# ------------------------------------------------------------------------------------ P4/P5
def p4(args):
    (F, C), _ = sizes(args)
    fl = FL(F, C)
    header(f"P4  facility location, F={F} C={C}: {fl.ncol} cols, {fl.nrow} rows, {fl.nnz} nnz")
    reps = 1 if args.quick else 3
    ns, how = build_rusolve(rs_expr, fl, reps)
    report("P4", f"rusolve expression API (build via {how})", ns, "ns/nnz", 300.0)
    P4["expr"] = ns
    if not args.compare:
        return
    for name, fn in (("highspy expression API", hs_expr), ("pyoptinterface", poi_loop)):
        try:
            t = timed(lambda fn=fn: fn(fl), 1)
            report("P4", f"  {name}", t / fl.nnz / NS, "ns/nnz")
        except ImportError:
            print(f"P4  {name} missing, row skipped")
    rows = dict_rows(fl)
    ns, _ = build_rusolve(lambda f: rs_dict(f, rows), fl, reps)
    report("P4", "  rusolve dict rows (dicts prebuilt)", ns, "ns/nnz")
    ns, _ = build_rusolve(rs_csr, fl, reps)
    report("P4", "  rusolve numpy CSR + add_matrix (assembly incl.)", ns, "ns/nnz")
    try:
        t = timed(lambda: hs_csr(fl), reps)
        report("P4", "  highspy numpy CSR addRows (assembly incl.)", t / fl.nnz / NS, "ns/nnz")
    except ImportError:
        print("P4  highspy missing, row skipped")


P4 = {}


def p5(args):
    (F, C), (F2, C2) = sizes(args)
    small, big = FL(F, C), FL(F2, C2)
    header(f"P5  scaling: ns/nnz at {big.nnz} nnz within 1.5x of the value at {small.nnz} nnz")
    reps = 1 if args.quick else 2
    prebuilt = {}

    def dict_path(fl):
        if fl not in prebuilt:
            prebuilt[fl] = dict_rows(fl)
        return rs_dict(fl, prebuilt[fl])

    paths = (("expression API", rs_expr), ("dict rows (dicts prebuilt)", dict_path),
             ("numpy CSR + add_matrix", rs_csr))
    for name, fn in paths:
        prebuilt.clear()
        a, _ = build_rusolve(fn, small, reps)
        gc.collect()
        b, _ = build_rusolve(fn, big, reps)
        prebuilt.clear()
        report("P5", f"{name}: {a:.0f} -> {b:.0f} ns/nnz, ratio", b / a, "x", 1.5)


# ------------------------------------------------------------------------------------ P6
def p6(args):
    n = 200_000 if args.quick else 1_000_000
    header(f"P6  sol.x for {n} variables")
    m = rs.Model()
    m.add_vars(n, ub=1, cost=-1.0)
    sol = m.solve()
    assert sol.ok, sol.status
    reps = 20 if args.quick else 50
    ts = []
    for _ in range(reps):
        s = m.solve()  # fresh Solution each time: first access includes the conversion
        t = time.perf_counter_ns()
        x = s.x
        ts.append(time.perf_counter_ns() - t)
    assert len(x) == n
    ms = statistics.median(ts) / 1e6
    report("P6", f"sol.x, first access, {n} vars (numpy: {isinstance(x, np.ndarray)})", ms, "ms",
           2.0 * n / 1_000_000)
    t = time.perf_counter_ns()
    _ = sol.x
    _ = sol.x
    report("P6", "  sol.x, two accesses of the same Solution", (time.perf_counter_ns() - t) / 1e6,
           "ms")


# ----------------------------------------------------------------------------------- main
ALL = {"p1": p1, "p2": p2, "p3": p3, "p4": p4, "p5": p5, "p6": p6}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("which", nargs="*", metavar="Pn")
    ap.add_argument("--quick", action="store_true", help="short run, small sizes")
    ap.add_argument("--no-compare", dest="compare", action="store_false",
                    help="skip highspy / pyoptinterface / scipy comparison rows")
    args = ap.parse_args()
    which = [w.lower() for w in args.which] or list(ALL)
    bad = [w for w in which if w not in ALL]
    if bad:
        ap.error(f"unknown benchmark {bad}; choose from {list(ALL)}")
    print(f"rusolve {getattr(rs, '__version__', '?')}  python {sys.version.split()[0]}  "
          f"{'quick' if args.quick else 'full'}")
    for w in which:
        ALL[w](args)
    fails = [r for r in RESULTS if not r[5]]
    print(f"\n{len(RESULTS) - len(fails)} PASS, {len(fails)} FAIL")
    for pid, what, value, unit, target, _ in fails:
        print(f"  FAIL {pid} {what}: {value:.2f} {unit} > {target:g}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
