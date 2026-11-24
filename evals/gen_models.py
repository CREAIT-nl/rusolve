#!/usr/bin/env python3
"""Writes the eval model files (deterministic: a fixed seed per model).

Run from anywhere: `python3 evals/gen_models.py`. The committed MPS files are the
output of this script; the graders never trust it, they recompute ground truth with
rusolve itself.
"""
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


def write_mps(path, name, cols, rows, rhs, bounds=None, int_cols=()):
    """cols: {col: {row_or_obj: coef}}; rows: {row: 'L'|'G'|'E'}; the objective row is 'obj'."""
    out = [f"NAME {name}", "ROWS", " N obj"]
    out += [f" {t} {r}" for r, t in rows.items()]
    out.append("COLUMNS")
    in_int = False
    for c, coefs in cols.items():
        if c in int_cols and not in_int:
            out.append("    MARKER  'MARKER'  'INTORG'")
            in_int = True
        if c not in int_cols and in_int:
            out.append("    MARKER  'MARKER'  'INTEND'")
            in_int = False
        for r, v in coefs.items():
            out.append(f"    {c} {r} {v:g}")
    if in_int:
        out.append("    MARKER  'MARKER'  'INTEND'")
    out.append("RHS")
    out += [f"    rhs {r} {v:g}" for r, v in rhs.items()]
    if bounds:
        out.append("BOUNDS")
        out += [f" {t} bnd {c} {v:g}" for c, t, v in bounds]
    out.append("ENDATA")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def multi_knapsack(path, name, n, m, seed):
    rnd = random.Random(seed)
    w = [[rnd.randint(5, 40) for _ in range(n)] for _ in range(m)]
    val = [sum(w[i][j] for i in range(m)) // m + rnd.randint(0, 9) for j in range(n)]
    cols = {f"x{j}": {"obj": -val[j], **{f"k{i}": w[i][j] for i in range(m)}} for j in range(n)}
    rows = {f"k{i}": "L" for i in range(m)}
    rhs = {f"k{i}": sum(w[i]) // 2 for i in range(m)}
    bounds = [(f"x{j}", "BV", 1) for j in range(n)]
    write_mps(path, name, cols, rows, rhs, bounds, int_cols=set(cols))


def single_knapsack(path, name, weights, profits, cap):
    cols = {f"x{j}": {"obj": -profits[j], "cap": weights[j]} for j in range(len(weights))}
    bounds = [(c, "BV", 1) for c in cols]
    write_mps(path, name, cols, {"cap": "L"}, {"cap": cap}, bounds, int_cols=set(cols))


def infeasible_plan(path):
    """Row `dem1` needs a+b >= 30 but a <= 15 and b <= 10: infeasible. The other rows are decoys."""
    cols = {
        "a": {"obj": 2, "dem1": 1, "cap1": 1},
        "b": {"obj": 3, "dem1": 1, "cap2": 1},
        "c": {"obj": 1, "dem2": 1, "cap1": 1},
        "d": {"obj": 4, "dem2": 1, "cap2": 1},
    }
    rows = {"dem1": "G", "dem2": "G", "cap1": "L", "cap2": "L"}
    rhs = {"dem1": 30, "dem2": 25, "cap1": 40, "cap2": 40}
    bounds = [("a", "UP", 15), ("b", "UP", 10), ("c", "UP", 20), ("d", "UP", 20)]
    write_mps(path, "INFPLAN", cols, rows, rhs, bounds)


def verify_model(path):
    cols = {
        "x": {"obj": -3, "cap": 2, "dem": 1, "bal": 1},
        "y": {"obj": -2, "cap": 3, "dem": 1, "bal": -1},
        "z": {"obj": -4, "cap": 4, "dem": 1},
        "w": {"obj": -1, "cap": 1},
    }
    rows = {"cap": "L", "dem": "G", "bal": "E"}
    rhs = {"cap": 20, "dem": 6, "bal": 0}
    bounds = [("x", "UP", 6), ("y", "UP", 6), ("z", "UP", 6), ("w", "UP", 8)]
    write_mps(path, "VERIFY", cols, rows, rhs, bounds, int_cols={"x", "y", "z"})


# ---------------------------------------------------------------- problem-class tasks
# One small model per playbook in docs/skills (LP file format: it carries indicators, semi-continuous
# columns, piecewise-linear constraints and a quadratic objective, which MPS here cannot). The
# optimum is never stored: graders.py recomputes it with rusolve.

def write_lp(path, sense, objective, rows, bounds=(), binaries=(), generals=(), semis=(), gcs=()):
    """objective: text after `obj:`; rows: [(name, text)]; bounds: [text]; gcs: general-constraint lines."""
    out = [sense, f" obj: {objective}", "Subject To"] + [f" {n}: {t}" for n, t in rows]
    if bounds:
        out += ["Bounds"] + [f" {b}" for b in bounds]
    if generals:
        out += ["General", " " + " ".join(generals)]
    if binaries:
        out += ["Binaries", " " + " ".join(binaries)]
    if semis:
        out += ["Semi-continuous", " " + " ".join(semis)]
    if gcs:
        out += ["General Constraints"] + [f" {g}" for g in gcs]
    out.append("End")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def _sum(terms):
    return " + ".join(terms)


def class_flow(path):
    """Transportation: 3 plants, 4 markets, balanced (100 units)."""
    supply, demand = [30, 40, 30], [20, 25, 30, 25]
    cost = [[8, 6, 10, 9], [9, 12, 13, 7], [14, 9, 16, 5]]
    x = lambda i, j: f"ship_{i}_{j}"
    obj = _sum(f"{cost[i][j]} {x(i, j)}" for i in range(3) for j in range(4))
    rows = [(f"supply_{i}", _sum(x(i, j) for j in range(4)) + f" = {supply[i]}") for i in range(3)]
    rows += [(f"demand_{j}", _sum(x(i, j) for i in range(3)) + f" = {demand[j]}") for j in range(4)]
    write_lp(path, "Minimize", obj, rows)


def class_cover(path):
    """Set cover: 12 elements, 9 sets with costs; every element in at least one set."""
    sets = {"S0": ([0, 1, 2, 3], 5), "S1": ([2, 3, 4, 5], 4), "S2": ([4, 5, 6, 7], 5), "S3": ([6, 7, 8, 9], 4),
            "S4": ([8, 9, 10, 11], 5), "S5": ([0, 4, 8], 3), "S6": ([1, 5, 9], 3), "S7": ([2, 6, 10, 11], 4),
            "S8": ([3, 7, 11, 0], 4)}
    obj = _sum(f"{c} {s}" for s, (_, c) in sets.items())
    rows = [(f"cover_{e}", _sum(s for s, (els, _) in sets.items() if e in els) + " >= 1") for e in range(12)]
    write_lp(path, "Minimize", obj, rows, binaries=list(sets))


def class_scheduling(path):
    """One machine, release dates, minimise makespan; disjunctions with a deliberately weak big-M
    (1000 against a horizon of 40): the classic mistake the playbook teaches lint to find."""
    p, r = [6, 4, 7, 3, 5], [0, 2, 1, 8, 4]
    n, H, M = 5, 40, 1000
    rows = []
    for i in range(n):
        rows.append((f"mk_{i}", f"cmax - s{i} >= {p[i]}"))
    for i in range(n):
        for j in range(i + 1, n):
            rows.append((f"ij_{i}_{j}", f"s{i} - s{j} + {M} y_{i}_{j} <= {M - p[i]}"))
            rows.append((f"ji_{i}_{j}", f"s{j} - s{i} - {M} y_{i}_{j} <= {-p[j]}"))
    bounds = [f"{r[i]} <= s{i} <= {H}" for i in range(n)] + [f"0 <= cmax <= {H + 20}"]
    write_lp(path, "Minimize", "cmax", rows, bounds, binaries=[f"y_{i}_{j}" for i in range(n) for j in range(i + 1, n)])


def class_lotsizing(path):
    """Capacitated single-item lot sizing, 8 periods, setup 90, linking x <= 5000 y (weak)."""
    d = [20, 40, 0, 35, 50, 15, 45, 30]
    T, cap, M = len(d), 80, 5000
    obj = _sum([f"2 x{t} + 1 i{t} + 90 y{t}" for t in range(T)])
    rows = []
    for t in range(T):
        prev = f"i{t - 1} + " if t else ""
        rows.append((f"bal_{t}", f"{prev}x{t} - i{t} = {d[t]}"))
        rows.append((f"link_{t}", f"x{t} - {M} y{t} <= 0"))
    bounds = [f"0 <= x{t} <= {cap}" for t in range(T)] + [f"0 <= i{t} <= 200" for t in range(T)]
    write_lp(path, "Minimize", obj, rows, bounds, binaries=[f"y{t}" for t in range(T)])


def class_facility(path):
    """Capacitated facility location, 4 sites, 6 customers, split demand."""
    fixed, cap = [60, 70, 50, 65], [45, 50, 40, 50]
    dem = [10, 12, 8, 15, 9, 11]
    c = [[4, 6, 9, 5, 8, 7], [7, 3, 5, 6, 4, 9], [6, 8, 3, 7, 5, 4], [5, 7, 6, 3, 9, 6]]
    obj = _sum([f"{fixed[j]} open_{j}" for j in range(4)] + [f"{dem[i] * c[j][i]} x_{j}_{i}" for j in range(4) for i in range(6)])
    rows = [(f"serve_{i}", _sum(f"x_{j}_{i}" for j in range(4)) + " = 1") for i in range(6)]
    rows += [(f"cap_{j}", _sum(f"{dem[i]} x_{j}_{i}" for i in range(6)) + f" - {cap[j]} open_{j} <= 0") for j in range(4)]
    rows += [(f"link_{j}_{i}", f"x_{j}_{i} - open_{j} <= 0") for j in range(4) for i in range(6)]
    bounds = [f"0 <= x_{j}_{i} <= 1" for j in range(4) for i in range(6)]
    write_lp(path, "Minimize", obj, rows, bounds, binaries=[f"open_{j}" for j in range(4)])


def class_routing(path):
    """Symmetric TSP on 6 nodes (node 0 is the depot), MTZ subtour elimination."""
    pts = [(0, 0), (4, 1), (7, 5), (3, 8), (-2, 5), (-3, 1)]
    n = len(pts)
    dist = lambda a, b: round(((pts[a][0] - pts[b][0]) ** 2 + (pts[a][1] - pts[b][1]) ** 2) ** 0.5 * 10)
    arcs = [(i, j) for i in range(n) for j in range(n) if i != j]
    obj = _sum(f"{dist(i, j)} x_{i}_{j}" for i, j in arcs)
    rows = [(f"out_{i}", _sum(f"x_{i}_{j}" for j in range(n) if j != i) + " = 1") for i in range(n)]
    rows += [(f"in_{j}", _sum(f"x_{i}_{j}" for i in range(n) if i != j) + " = 1") for j in range(n)]
    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                rows.append((f"mtz_{i}_{j}", f"u{i} - u{j} + {n - 1} x_{i}_{j} <= {n - 2}"))
    bounds = [f"1 <= u{i} <= {n - 1}" for i in range(1, n)]
    write_lp(path, "Minimize", obj, rows, bounds, binaries=[f"x_{i}_{j}" for i, j in arcs])


def class_pwl_semi(path):
    """Two plants must supply 60 units. Plant a: semi-continuous output (off, or 15..40) and a
    non-convex piecewise-linear cost (cheap first block, a surcharge, then cheap again); plant b:
    plain linear cost 3.1 per unit up to 50."""
    rows = [("demand", "a + b >= 60")]
    bounds = ["15 <= a <= 40", "0 <= b <= 50", "0 <= ca <= 400"]
    gcs = ["pa: ca = PWL ( a ) : (0, 0) (15, 30) (25, 90) (40, 110)"]
    write_lp(path, "Minimize", "ca + 3.1 b", rows, bounds, semis=["a"], gcs=gcs)


def class_qp(path):
    """Mean-variance portfolio: 4 assets, full investment, expected return >= 0.08, no shorting."""
    mu = [0.05, 0.08, 0.11, 0.14]
    cov = [[0.04, 0.006, 0.002, 0.0], [0.006, 0.09, 0.012, 0.01], [0.002, 0.012, 0.16, 0.03], [0.0, 0.01, 0.03, 0.25]]
    q = []
    for i in range(4):
        for j in range(i, 4):
            q.append(f"{cov[i][j] * 2:g} w{i} * w{j}" if i != j else f"{cov[i][i]:g} w{i} ^ 2")
    obj = "[ " + " + ".join(q) + " ] / 2"
    rows = [("budget", _sum(f"w{i}" for i in range(4)) + " = 1"), ("ret", _sum(f"{mu[i]} w{i}" for i in range(4)) + " >= 0.08")]
    write_lp(path, "Minimize", obj, rows, [f"0 <= w{i} <= 1" for i in range(4)])


def class_conic(path):
    """Revenue maximisation under an ellipsoidal risk budget: a convex quadratic row (a second-order
    cone constraint written as a quadratic row) and a linear budget; no integers."""
    rows = [("risk", "[ 0.04 a ^ 2 + 0.09 b ^ 2 + 0.16 c ^ 2 + 0.02 a * b ] <= 1"), ("budget", "a + b + c <= 10")]
    write_lp(path, "Maximize", "4 a + 3 b + 2 c", rows, ["0 <= a <= 10", "0 <= b <= 10", "0 <= c <= 10"])


def class_miqp(path):
    """Cardinality-constrained portfolio (convex MIQP): the class_qp assets, risk aversion 10, at most two
    assets held, each at least 10 % when held; the objective is -mu'w + 5 w'cov w."""
    mu = [0.05, 0.08, 0.11, 0.14]
    cov = [[0.04, 0.006, 0.002, 0.0], [0.006, 0.09, 0.012, 0.01], [0.002, 0.012, 0.16, 0.03], [0.0, 0.01, 0.03, 0.25]]
    q = []
    for i in range(4):
        for j in range(i, 4):
            q.append(f"{cov[i][j] * 20:g} w{i} * w{j}" if i != j else f"{cov[i][i] * 10:g} w{i} ^ 2")
    obj = " ".join(f"- {mu[i]} w{i}" for i in range(4)) + " + [ " + " + ".join(q) + " ] / 2"
    rows = [("budget", _sum(f"w{i}" for i in range(4)) + " = 1"), ("cardinality", _sum(f"y{i}" for i in range(4)) + " <= 2")]
    rows += [(f"hi_{i}", f"w{i} - y{i} <= 0") for i in range(4)]
    rows += [(f"lo_{i}", f"w{i} - 0.1 y{i} >= 0") for i in range(4)]
    write_lp(path, "Minimize", obj, rows, [f"0 <= w{i} <= 1" for i in range(4)], binaries=[f"y{i}" for i in range(4)])


def class_misocp(path):
    """Mixed-integer convex quadratic rows: whole batches a, b and a fractional amount c, an energy
    row a^2 + 2 b^2 + c^2 <= 40 and a weight row; maximise 5a + 4b + 3c."""
    rows = [("energy", "[ a ^ 2 + 2 b ^ 2 + c ^ 2 ] <= 40"), ("weight", "a + b + c <= 8")]
    write_lp(path, "Maximize", "5 a + 4 b + 3 c", rows, ["0 <= a <= 6", "0 <= b <= 6", "0 <= c <= 10"], generals=["a", "b"])


def supplier_sourcing_reference(path):
    """Reference model of the word problem `supplier_sourcing` (general constraints), written by hand
    with the names the prompt dictates."""
    rows = [("demand", "qA + qB + qC >= 100")]
    gcs = ["g_largest: largest = MAX ( qA , qB , qC )"]
    out = ["Minimize", " obj: 9 qA + 7 qB + 8 qC + 120 useB + 0.2 largest", "Subject To",
           " demand: qA + qB + qC >= 100", " b_off: useB = 0 -> qB <= 0", " b_min: useB = 1 -> qB >= 25",
           "Bounds", " 0 <= qA <= 60", " 0 <= qB <= 70", " 0 <= qC <= 50", " 0 <= largest <= 70",
           "Binaries", " useB", "SOS", " a_or_c: S1 :: qA : 1 qC : 2", "General Constraints"] + [f" {g}" for g in gcs] + ["End"]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    t = os.path.join(HERE, "tasks")
    infeasible_plan(f"{t}/debug_infeasible/model.mps")
    verify_model(f"{t}/verify_formulation/model.mps")
    single_knapsack(f"{t}/add_cuts/model.mps", "CUTKNAP",
                    [34, 33, 32, 31, 30, 29, 28, 27, 26, 25], [36, 34, 34, 32, 31, 31, 29, 28, 27, 27], 100)
    # solve_report: 12-item two-constraint knapsack.
    multi_knapsack(f"{t}/solve_report/model.mps", "SOLVEREP", 12, 2, 11)
    # tune_family: three instances of one family (same generator, different seeds).
    for s in (1, 2, 3):
        multi_knapsack(f"{t}/tune_family/inst{s}.mps", f"FAM{s}", 40, 5, 100 + s)
    # extend_wasm: 30 items, 4 constraints, seed 7. The search with the HiGHS child rule needs 912
    # nodes here (measured, `rusolve mip --format json --seed 0 --opt mip_child_selection=-1`), so a
    # node-selection rule has a real tree to improve; the evals README records the numbers.
    multi_knapsack(f"{t}/extend_wasm/model.mps", "WASMEXT", 30, 4, 7)
    # Problem-class tasks (docs/skills): one model per playbook.
    for name, fn in [("flow", class_flow), ("cover", class_cover), ("scheduling", class_scheduling),
                     ("lotsizing", class_lotsizing), ("facility", class_facility), ("routing", class_routing),
                     ("pwl_semi", class_pwl_semi), ("qp", class_qp), ("conic", class_conic),
                     ("miqp", class_miqp), ("misocp", class_misocp)]:
        fn(f"{t}/class_{name}/model.lp")
    supplier_sourcing_reference(f"{t}/supplier_sourcing/selfcheck/reference.lp")


if __name__ == "__main__":
    main()
