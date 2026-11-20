"""Deterministic instance generators for differential testing against HiGHS and SCIP.

Every generator returns the generic intermediate model defined below and is
seeded (default seed 42) so runs are bit-reproducible. Default sizes are small
(at most 1000 variables and 1000 rows), so every instance solves quickly.
Used by tools/diff/mip_diff.py.

Provenance of each family: standard benchmark literature (degenerate LPs,
set-covering, transshipment, knapsack-cover, Klee-Minty), cited per family;
the mixed-structure families follow the row classes of Gleixner et al.,
"MIPLIB 2017", Math. Prog. Comp. 13 (2021).
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

INF = float("inf")

VTYPE_CONT = "continuous"
VTYPE_INT = "integer"
VTYPE_BIN = "binary"


@dataclass(frozen=True)
class Var:
    lb: float = 0.0
    ub: float = INF
    vtype: str = VTYPE_CONT
    name: str = ""


@dataclass(frozen=True)
class Row:
    """Linear row: sum(coeff * x[idx]) sense rhs."""

    coeffs: Tuple[Tuple[int, float], ...]  # (var_index, coefficient) pairs
    sense: str  # "<=", ">=", "=="
    rhs: float
    name: str = ""


@dataclass
class Model:
    """Solver-agnostic intermediate model."""

    name: str
    vars: List[Var]
    rows: List[Row]
    obj: List[Tuple[int, float]]  # (var_index, coefficient)
    obj_sense: str = "minimize"  # "minimize" | "maximize"
    meta: Dict[str, object] = field(default_factory=dict)

    @property
    def ncols(self) -> int:
        return len(self.vars)

    @property
    def nrows(self) -> int:
        return len(self.rows)

    @property
    def is_mip(self) -> bool:
        return any(v.vtype != VTYPE_CONT for v in self.vars)


# ---------------------------------------------------------------------------


def validate_planted(mdl: "Model", eps: float = 1e-7) -> "str | None":
    """Check a generator's planted witness solution against bounds/rows/objective.

    Generators that can guarantee feasibility attach meta['planted'] (a
    var-index -> value mapping or full-length list); families with a known
    optimum also attach meta['planted_obj'. Returns None when the witness is
    valid (or none was planted), else a human-readable mismatch description.
    Catches generator bugs at generation time instead of relying on solvers
    agreeing on 'infeasible'.
    """
    planted = mdl.meta.get("planted")
    if planted is None:
        return None
    x = [0.0] * mdl.ncols
    if isinstance(planted, dict):
        for j, v in planted.items():
            x[int(j)] = float(v)
    else:
        for j, v in enumerate(planted):
            x[j] = float(v)
    for j, v in enumerate(mdl.vars):
        if x[j] < v.lb - eps or x[j] > v.ub + eps:
            return f"var {mdl.vars[j].name or j}={x[j]} violates bounds [{v.lb}, {v.ub}]"
    for r in mdl.rows:
        lhs = sum(c * x[int(i)] for i, c in r.coeffs)
        if r.sense == "<=" and lhs > r.rhs + eps:
            return f"row {r.name}: {lhs} > rhs {r.rhs}"
        if r.sense == ">=" and lhs < r.rhs - eps:
            return f"row {r.name}: {lhs} < rhs {r.rhs}"
        if r.sense == "==" and abs(lhs - r.rhs) > eps:
            return f"row {r.name}: {lhs} != rhs {r.rhs}"
    pobj = mdl.meta.get("planted_obj")
    if pobj is not None:
        obj = sum(c * x[int(j)] for j, c in mdl.obj)
        if abs(obj - float(pobj)) > eps * max(1.0, abs(float(pobj))):
            return f"planted obj {obj} != predicted {pobj}"
    return None


def degenerate_ladder(n: int = 50, seed: int = 42) -> Model:
    """LP family with massively degenerate vertices and tied objective costs.

    Structure (2n vars, n rows):
        x_i + y_i = 1        i = 1..n
        x_{i+1} - y_i = 0    i = 1..n-1
        0 <= x_i, y_i <= 1
    The feasible set is a 1-D segment (x alternates t, 1-t); both endpoints are
    vertices where ALL 2n bounds are active while the basis has size n ->
    maximal degeneracy. Objective costs are drawn from {1, 2} to create massive
    reduced-cost ties along the way. Probes dual-ratio ties / perturbation.
    """
    rng = random.Random(seed)
    vs: List[Var] = []
    rows: List[Row] = []
    for i in range(n):
        vs.append(Var(name=f"x{i}"))
        vs.append(Var(name=f"y{i}"))
    for i in range(n):
        rows.append(Row(((2 * i, 1.0), (2 * i + 1, 1.0)), "==", 1.0, name=f"cover{i}"))
    for i in range(n - 1):
        rows.append(Row(((2 * (i + 1), 1.0), (2 * i + 1, -1.0)), "==", 0.0, name=f"link{i}"))
    obj = [(2 * i, float(rng.choice((1, 2)))) for i in range(n)]
    # Feasibility witness: x_i = y_i = 0.5 satisfies every row and bound.
    return Model(f"degenerate_ladder_{n}", vs, rows, obj, "minimize",
                 meta={"family": "degenerate_ladder", "n": n, "seed": seed,
                       "planted": [0.5] * (2 * n)})


def set_cover(m: int = 60, n: int = 120, density: float = 0.1, seed: int = 42) -> Model:
    """Random set-covering MIP: min c'x s.t. Ax >= 1, x binary.

    Each column covers each row independently with prob `density`; rows with no
    cover are repaired. Costs are small integers -> many ties.
    """
    rng = random.Random(seed)
    vs = [Var(0, 1, VTYPE_BIN, name=f"x{j}") for j in range(n)]
    cols: List[List[int]] = [[] for _ in range(n)]
    for j in range(n):
        for i in range(m):
            if rng.random() < density:
                cols[j].append(i)
    covered = [any(i in c for c in cols) for i in range(m)]
    for i in range(m):
        if not covered[i]:
            j = rng.randrange(n)
            cols[j].append(i)
    rows = []
    for i in range(m):
        coeffs = tuple((j, 1.0) for j in range(n) if i in cols[j])
        rows.append(Row(coeffs, ">=", 1.0, name=f"cover{i}"))
    obj = [(j, float(rng.randint(1, 100))) for j in range(n)]
    # Feasibility witness: all-ones covers every row (empty rows were repaired).
    return Model(f"set_cover_{m}x{n}", vs, rows, obj, "minimize",
                 meta={"family": "set_cover", "m": m, "n": n, "density": density,
                       "seed": seed, "planted": [1] * n})


def network_flow(nodes: int = 40, arcs: int = 120, seed: int = 42) -> Model:
    """Min-cost transshipment LP with guaranteed feasible pure-network structure.

    Feasibility by construction: a random bidirectional spanning tree is built,
    integer arc flows are planted on it, and node balances are derived from the
    planted flows; extra random arcs are then added.
        min  sum c_a f_a
        s.t. net outflow at node i == b_i (all rows '==', sums to zero)
             0 <= f_a <= cap_a
    """
    rng = random.Random(seed)
    tree: List[Tuple[int, int]] = []
    for i in range(1, nodes):
        tree.append((i, rng.randrange(i)))
    arc_set = set()
    for u, v in tree:
        arc_set.add((u, v))
        arc_set.add((v, u))
    arcs_list = []
    caps = {}
    for u, v in tree:
        for a in ((u, v), (v, u)):
            caps[a] = rng.randint(3, 15)
            arcs_list.append(a)
    guard = 0
    while len(arcs_list) < arcs and guard < 10 * arcs:
        guard += 1
        u, v = rng.randrange(nodes), rng.randrange(nodes)
        if u != v and (u, v) not in arc_set:
            arc_set.add((u, v))
            caps[(u, v)] = rng.randint(3, 15)
            arcs_list.append((u, v))
    # Plant integer flows on the tree (both directions available, cap >= 3).
    bal = [0] * nodes
    planted_flows: Dict[Tuple[int, int], int] = {}
    for u, v in tree:
        # f was once bounded only by caps[(u,v)]//2 but the coin
        # below may plant it on arc (v,u), whose cap is drawn independently ->
        # planted flow could exceed cap[(v,u)] and the LP became infeasible
        # (seen at nodes=200/1000/2000, seed 42). Bound by BOTH caps.
        # (The infeasibility went unnoticed until solvers *agreed* on it —
        # validate_planted() now catches this class of bug at generation.)
        f = rng.randint(0, min(caps[(u, v)], caps[(v, u)]) // 2)
        if rng.random() < 0.5:
            bal[u] += f
            bal[v] -= f
            planted_flows[(u, v)] = planted_flows.get((u, v), 0) + f
        else:
            bal[v] += f
            bal[u] -= f
            planted_flows[(v, u)] = planted_flows.get((v, u), 0) + f
    vs = [Var(0, float(caps[a]), name=f"f{a[0]}_{a[1]}") for a in arcs_list]
    idx = {a: k for k, a in enumerate(arcs_list)}
    rows = []
    for i in range(nodes):
        coeffs = []
        for k, (u, v) in enumerate(arcs_list):
            if u == i:
                coeffs.append((k, 1.0))
            elif v == i:
                coeffs.append((k, -1.0))
        rows.append(Row(tuple(coeffs), "==", float(bal[i]), name=f"bal{i}"))
    obj = [(idx[a], float(rng.randint(1, 100))) for a in arcs_list]
    # Planted flow: tree arcs carry the planted amounts (within caps by the
    # min(both caps) bound above), all other arcs 0; balances are derived from
    # exactly these flows, so the witness satisfies every row at construction.
    planted = {idx[a]: float(val) for a, val in planted_flows.items()}
    return Model(f"network_flow_{nodes}n", vs, rows, obj, "minimize",
                 meta={"family": "network_flow", "nodes": nodes, "arcs": len(arcs_list),
                       "seed": seed, "total_supply": sum(b for b in bal if b > 0),
                       "planted": planted})


def knapsack_cover(n: int = 40, seed: int = 42) -> Model:
    """Hard knapsack-cover: min p'x s.t. w'x >= W, x binary (LP-relaxable).

    Inverse-strongly-correlated profits (p_i = W_MAX - w_i) make the LP
    relaxation weak -> many cover cuts needed. Classic hard-MIP family for
    exercising cut selection. Set meta['relax']=True (via --kw) for the LP version.
    """
    rng = random.Random(seed)
    w = [rng.randint(1, 100) for _ in range(n)]
    W = sum(w) // 2
    p = [101 - wi for wi in w]
    vs = [Var(0, 1, VTYPE_BIN, name=f"x{i}") for i in range(n)]
    rows = [Row(tuple((i, float(wi)) for i, wi in enumerate(w)), ">=", float(W), name="knap")]
    obj = [(i, float(p[i])) for i in range(n)]
    # Feasibility witness: all-ones (sum w = 2W >= W).
    return Model(f"knapsack_cover_{n}", vs, rows, obj, "minimize",
                 meta={"family": "knapsack_cover", "n": n, "seed": seed,
                       "capacity": W, "planted": [1] * n})


def corner_lp(n: int = 8, seed: int = 42) -> Model:
    """Klee-Minty cube: LP whose unique optimum is the all-tight far vertex.

        max sum_{j} 10^(n-j) x_j
        s.t. x_1 <= 1
             2*sum_{j<i} x_j + x_i <= 10^(2(i-1))   i = 2..n
             x >= 0
    A naive Dantzig rule can walk exponentially many vertices to reach the
    all-tight vertex x_1=1, x_i = 10^(2(i-1)) - 2*sum_{j<i} x_j (verified for
    n=2 by hand; HiGHS and SCIP agree on the value for n=8). The exact
    vertex/objective are in meta['predicted_vertex'] / ['predicted_obj'].
    """
    vs = [Var(0.0, INF, name=f"x{i}") for i in range(n)]
    rows: List[Row] = [Row(((0, 1.0),), "<=", 1.0, name="k0")]
    for i in range(1, n):
        coeffs = tuple([(j, 2.0) for j in range(i)] + [(i, 1.0)])
        rows.append(Row(coeffs, "<=", float(10 ** (2 * i)), name=f"k{i}"))
    obj = [(j, float(10 ** (n - 1 - j))) for j in range(n)]
    vertex = [1.0]
    for i in range(1, n):
        vertex.append(10.0 ** (2 * i) - 2.0 * sum(vertex))
    predicted_obj = sum((10.0 ** (n - 1 - j)) * vertex[j] for j in range(n))
    return Model(f"corner_lp_{n}", vs, rows, obj, "maximize",
                 meta={"family": "corner_lp", "n": n, "seed": seed,
                       "predicted_vertex": vertex,
                       "predicted_obj": predicted_obj,
                       "planted": vertex, "planted_obj": predicted_obj})


def transport_degen(m: int = 8, seed: int = 42) -> Model:
    """Presolve-RESISTANT degenerate LP: square transportation problem.

    Motivation: the degenerate_ladder's 2-nonzero equality chain is folded by
    a presolve aggregator, so the simplex never sees the degeneracy. Here every
    row has m >= 3 nonzeros, no row is singleton, every column lives in exactly
    two coupled equality rows (supply + demand) -> the classic aggregator has
    nothing to substitute, so the solver must run its LP algorithm.

    Degeneracy machinery:
      * integer margins s_i / d_j built from small planted flows x0_ij in
        {0,1,2} -> many equal margins -> degenerate transportation vertices
        (a BFS is degenerate whenever fewer than 2m-1 basic cells are positive);
      * costs drawn from {1,2,3} -> massive reduced-cost ties for pricing;
      * upper bounds u_ij = x0_ij + 1 keep planted flow strictly interior, so
        the planted witness proves feasibility (validated at generation time).
    Feasibility by construction: margins ARE the row/column sums of x0.
    """
    rng = random.Random(seed)
    x0 = [[rng.randint(0, 2) for _ in range(m)] for _ in range(m)]
    # guarantee every margin >= 1 so no row collapses to all-zero (which
    # presolve would use to fix whole rows of variables to 0)
    for i in range(m):
        if sum(x0[i]) == 0:
            x0[i][rng.randrange(m)] = 1
    for j in range(m):
        if sum(x0[i][j] for i in range(m)) == 0:
            x0[rng.randrange(m)][j] = 1
    supply = [sum(x0[i]) for i in range(m)]
    demand = [sum(x0[i][j] for i in range(m)) for j in range(m)]
    vs, obj, planted = [], [], {}
    for i in range(m):
        for j in range(m):
            k = i * m + j
            vs.append(Var(0.0, float(x0[i][j] + 1), name=f"x{i}_{j}"))
            obj.append((k, float(rng.choice((1, 2, 3)))))
            planted[k] = float(x0[i][j])
    rows = []
    for i in range(m):
        rows.append(Row(tuple((i * m + j, 1.0) for j in range(m)),
                        "==", float(supply[i]), name=f"sup{i}"))
    for j in range(m):
        rows.append(Row(tuple((i * m + j, 1.0) for i in range(m)),
                        "==", float(demand[j]), name=f"dem{j}"))
    return Model(f"transport_degen_{m}", vs, rows, obj, "minimize",
                 meta={"family": "transport_degen", "m": m, "seed": seed,
                       "planted": planted,
                       "supply_hist": {s: supply.count(s) for s in set(supply)},
                       "demand_hist": {d: demand.count(d) for d in set(demand)}})


# ---------------------------------------------------------------------------
# MIP families. Every family below:
#   * is deterministic per seed and has at most 1000 vars and 1000 rows at its
#     default size;
#   * plants a feasible witness in meta['planted'] where one is constructible
#     (checked by validate_planted);
#   * has a `small` size in MIP_SMALL (<= ~30 vars, <= ~25 rows) for a slow
#     exact B&B / the mip_diff correctness gate (tools/diff/mip_diff.py);
#   * infeasible=True variants (set_partition, gen_int_knapsack, gap) are
#     provably infeasible by a counting/parity argument while their LP
#     relaxation stays feasible -> the infeasibility must be proven by B&B.


def _bins(n: int, name: str = "x") -> List[Var]:
    return [Var(0, 1, VTYPE_BIN, name=f"{name}{j}") for j in range(n)]


def set_packing(m: int = 60, n: int = 120, density: float = 0.08, seed: int = 42) -> Model:
    """Random weighted set packing: max c'x s.t. Ax <= 1, x binary.

    Probes clique extraction / clique cuts and the clique table (every row is
    a clique); the LP relaxation is half-integral-prone. Formulation:
    Padberg 1973, "On the facial structure of set packing polyhedra",
    Math. Prog. 5. Witness: greedy packing in column order.
    """
    rng = random.Random(seed)
    cols: List[List[int]] = []
    for _ in range(n):
        c = [i for i in range(m) if rng.random() < density]
        if not c:
            c = [rng.randrange(m)]
        cols.append(c)
    rows = []
    for i in range(m):
        coeffs = tuple((j, 1.0) for j in range(n) if i in cols[j])
        rows.append(Row(coeffs, "<=", 1.0, name=f"pack{i}"))
    obj = [(j, float(rng.randint(1, 20))) for j in range(n)]
    used, planted = set(), [0] * n
    for j in range(n):
        if not used.intersection(cols[j]):
            used.update(cols[j])
            planted[j] = 1
    return Model(f"set_packing_{m}x{n}", _bins(n), rows, obj, "maximize",
                 meta={"family": "set_packing", "m": m, "n": n, "density": density,
                       "seed": seed, "planted": planted})


def set_partition(m: int = 40, n: int = 160, seed: int = 42, infeasible: bool = False) -> Model:
    """Set partitioning: min c'x s.t. Ax = 1, x binary.

    Probes clique cuts, the clique table, and the extreme dual degeneracy of
    partitioning LPs (airline crew-scheduling structure). Formulation:
    Balas & Padberg 1976, "Set partitioning: a survey", SIAM Review 18.
    Feasible mode: a random planted partition of the rows supplies columns
    that form an exact cover (meta['planted']); the rest are random subsets.
    infeasible=True: m is made odd and every column covers an EVEN number of
    rows, so 1'Ax = sum_j |S_j| x_j is even != m (parity proof), while
    x_j = 1/|S_j|-style fractional points keep the LP relaxation feasible.
    """
    rng = random.Random(seed)
    cols: List[List[int]] = []
    planted = None
    if infeasible:
        if m % 2 == 0:
            m += 1
        # Disjoint pairs + a closing pair cover every row exactly 2x at
        # x = 1/2 -> LP-feasible; every column has even cardinality.
        perm = list(range(m))
        rng.shuffle(perm)
        for k in range(0, m - 1, 2):
            cols.append(sorted((perm[k], perm[k + 1])))
        cols.append(sorted((perm[-1], perm[0])))
        for k in range(1, m - 1, 2):
            cols.append(sorted((perm[k], perm[k + 1])))
        while len(cols) < n:
            size = 2 * rng.randint(1, max(1, min(3, m // 2)))
            cols.append(sorted(rng.sample(range(m), size)))
    else:
        perm = list(range(m))
        rng.shuffle(perm)
        k = 0
        while k < m:
            size = rng.randint(1, 4)
            cols.append(sorted(perm[k:k + size]))
            k += size
        planted = [1] * len(cols)
        while len(cols) < n:
            cols.append(sorted(rng.sample(range(m), rng.randint(1, min(5, m)))))
        planted += [0] * (len(cols) - len(planted))
    # shuffle columns so the planted cover is not a prefix
    order = list(range(len(cols)))
    rng.shuffle(order)
    cols = [cols[k] for k in order]
    if planted is not None:
        planted = [planted[k] for k in order]
    n = len(cols)
    rows = []
    for i in range(m):
        rows.append(Row(tuple((j, 1.0) for j in range(n) if i in cols[j]), "==", 1.0,
                        name=f"part{i}"))
    obj = [(j, float(rng.randint(1, 10) * len(cols[j]) + rng.randint(0, 5))) for j in range(n)]
    meta = {"family": "set_partition", "m": m, "n": n, "seed": seed, "infeasible": infeasible}
    if planted is not None:
        meta["planted"] = planted
    tag = "_infeas" if infeasible else ""
    return Model(f"set_partition_{m}x{n}{tag}", _bins(n), rows, obj, "minimize", meta=meta)


def multi_knapsack(n: int = 100, m: int = 10, tightness: float = 0.5, seed: int = 42) -> Model:
    """Multidimensional 0-1 knapsack (MKP): max p'x s.t. Wx <= b, x binary.

    Probes lifted cover cuts, reduced-cost fixing and dive heuristics.
    Generator: Chu & Beasley 1998, "A genetic algorithm for the
    multidimensional knapsack problem", J. Heuristics 4 — w_ij ~ U[1,1000],
    b_i = tightness * sum_j w_ij, p_j = sum_i w_ij / m + 500 u_j (correlated).
    Witness: greedy by profit/weight ratio.
    """
    rng = random.Random(seed)
    w = [[rng.randint(1, 1000) for _ in range(n)] for _ in range(m)]
    b = [int(tightness * sum(w[i])) for i in range(m)]
    p = [int(sum(w[i][j] for i in range(m)) / m + 500 * rng.random()) for j in range(n)]
    rows = [Row(tuple((j, float(w[i][j])) for j in range(n)), "<=", float(b[i]), name=f"cap{i}")
            for i in range(m)]
    obj = [(j, float(p[j])) for j in range(n)]
    load, planted = [0] * m, [0] * n
    for j in sorted(range(n), key=lambda j: -p[j] / sum(w[i][j] for i in range(m))):
        if all(load[i] + w[i][j] <= b[i] for i in range(m)):
            planted[j] = 1
            for i in range(m):
                load[i] += w[i][j]
    return Model(f"multi_knapsack_{n}x{m}", _bins(n), rows, obj, "maximize",
                 meta={"family": "multi_knapsack", "n": n, "m": m, "tightness": tightness,
                       "seed": seed, "planted": planted})


def gen_int_knapsack(n: int = 60, m: int = 3, equality: bool = False, seed: int = 42,
                     infeasible: bool = False) -> Model:
    """Bounded general-integer knapsack: max p'x s.t. Wx <= b, 0 <= x_j <= u_j int.

    Probes MIR / c-MIR on general integers and branching on non-binaries
    (Marchand & Wolsey 2001, "Aggregation and mixed integer rounding to solve
    MIPs", Oper. Res. 49; bounded knapsack: Martello & Toth 1990, "Knapsack
    Problems", ch. 3). equality=True turns row 0 into w'x = b0 with b0 = w'x0
    for a planted integer x0 (the equality-knapsack lattice structure).
    infeasible=True: row 0 is an equality with ALL coefficients even and an
    ODD rhs -> no integer point (parity proof), yet the LP is feasible since
    the rhs lies strictly inside [0, w'u].
    """
    rng = random.Random(seed)
    u = [rng.randint(2, 10) for _ in range(n)]
    w = [[rng.randint(5, 60) for _ in range(n)] for _ in range(m)]
    if infeasible:
        w[0] = [2 * rng.randint(3, 30) for _ in range(n)]
    p = [int(sum(w[i][j] for i in range(m)) / m) + rng.randint(-4, 10) for j in range(n)]
    x0 = [rng.randint(0, u[j] // 2) for j in range(n)]
    b = [max(sum(w[i][j] * x0[j] for j in range(n)),
             int(0.3 * sum(w[i][j] * u[j] for j in range(n)))) for i in range(m)]
    eq0 = equality or infeasible
    if eq0:
        b[0] = sum(w[0][j] * x0[j] for j in range(n))
        if infeasible:
            b[0] += 1  # even LHS, odd rhs
    rows = []
    for i in range(m):
        sense = "==" if (i == 0 and eq0) else "<="
        rows.append(Row(tuple((j, float(w[i][j])) for j in range(n)), sense, float(b[i]),
                        name=f"knap{i}"))
    vs = [Var(0, float(u[j]), VTYPE_INT, name=f"x{j}") for j in range(n)]
    obj = [(j, float(p[j])) for j in range(n)]
    meta = {"family": "gen_int_knapsack", "n": n, "m": m, "equality": equality,
            "seed": seed, "infeasible": infeasible}
    if not infeasible:
        meta["planted"] = x0
    tag = "_infeas" if infeasible else ("_eq" if equality else "")
    return Model(f"gen_int_knapsack_{n}x{m}{tag}", vs, rows, obj, "maximize", meta=meta)


def _facility_data(rng, n_fac, n_cust):
    fx = [(rng.random(), rng.random()) for _ in range(n_fac)]
    cx = [(rng.random(), rng.random()) for _ in range(n_cust)]
    dist = [[((fx[i][0] - cx[j][0]) ** 2 + (fx[i][1] - cx[j][1]) ** 2) ** 0.5
             for j in range(n_cust)] for i in range(n_fac)]
    return dist


def cap_facility_location(n_fac: int = 12, n_cust: int = 40, ratio: float = 3.0,
                          strong: bool = True, seed: int = 42) -> Model:
    """Capacitated facility location (single-source relaxed, x_ij = fraction).

        min sum f_i y_i + sum c_ij d_j x_ij
        s.t. sum_i x_ij = 1                    (demand)
             sum_j d_j x_ij <= s_i y_i         (capacity, big-M-like coupling)
             x_ij <= y_i                       (strong linking, if strong)
             y binary, 0 <= x <= 1
    Probes flow-cover and implied-bound cuts, presolve strengthening. Instance
    generator after Cornuejols, Sridharan & Thizy 1991, "A comparison of
    heuristics and relaxations for the capacitated plant location problem",
    EJOR 50 (capacity ratio sum s / sum d = `ratio`). Witness: all facilities
    open, x_ij = s_i / sum s (each capacity row then holds iff sum d <= sum s).
    """
    rng = random.Random(seed)
    dist = _facility_data(rng, n_fac, n_cust)
    d = [rng.randint(5, 35) for _ in range(n_cust)]
    D = sum(d)
    s_raw = [rng.randint(10, 160) for _ in range(n_fac)]
    scale = ratio * D / sum(s_raw)
    s = [max(1, int(round(v * scale))) for v in s_raw]
    if sum(s) < D:
        s[0] += D - sum(s)
    f = [rng.randint(0, 90) + int(rng.randint(100, 110) * s[i] ** 0.5) for i in range(n_fac)]
    S = float(sum(s))
    vs = _bins(n_fac, "y")
    xi = lambda i, j: n_fac + i * n_cust + j  # noqa: E731
    for i in range(n_fac):
        for j in range(n_cust):
            vs.append(Var(0.0, 1.0, name=f"x{i}_{j}"))
    rows = []
    for j in range(n_cust):
        rows.append(Row(tuple((xi(i, j), 1.0) for i in range(n_fac)), "==", 1.0, name=f"dem{j}"))
    for i in range(n_fac):
        coeffs = tuple([(xi(i, j), float(d[j])) for j in range(n_cust)] + [(i, -float(s[i]))])
        rows.append(Row(coeffs, "<=", 0.0, name=f"cap{i}"))
    if strong:
        for i in range(n_fac):
            for j in range(n_cust):
                rows.append(Row(((xi(i, j), 1.0), (i, -1.0)), "<=", 0.0, name=f"link{i}_{j}"))
    obj = [(i, float(f[i])) for i in range(n_fac)]
    obj += [(xi(i, j), float(round(10 * dist[i][j] * d[j], 3)))
            for i in range(n_fac) for j in range(n_cust)]
    planted = [1.0] * n_fac + [s[i] / S for i in range(n_fac) for _ in range(n_cust)]
    return Model(f"cap_facility_location_{n_fac}x{n_cust}{'' if strong else '_weak'}", vs, rows,
                 obj, "minimize",
                 meta={"family": "cap_facility_location", "n_fac": n_fac, "n_cust": n_cust,
                       "ratio": ratio, "strong": strong, "seed": seed, "planted": planted})


def uncap_facility_location(n_fac: int = 20, n_cust: int = 40, form: str = "strong",
                            seed: int = 42) -> Model:
    """Uncapacitated facility location, weak vs strong linking.

        weak:   sum_j x_ij <= n_cust * y_i     (one aggregated row per facility)
        strong: x_ij <= y_i                    (disaggregated)
    The two forms have the same integer hull but very different LP bounds;
    probes presolve/cut strengthening of the weak form (implied-bound cuts
    recover the strong form). Cornuejols, Nemhauser & Wolsey 1990, "The
    uncapacitated facility location problem", in Discrete Location Theory.
    Witness: open facility 0 only, assign every customer to it.
    """
    if form not in ("weak", "strong"):
        raise ValueError("form must be 'weak' or 'strong'")
    rng = random.Random(seed)
    dist = _facility_data(rng, n_fac, n_cust)
    f = [rng.randint(20, 60) for _ in range(n_fac)]
    vs = _bins(n_fac, "y")
    xi = lambda i, j: n_fac + i * n_cust + j  # noqa: E731
    for i in range(n_fac):
        for j in range(n_cust):
            vs.append(Var(0.0, 1.0, name=f"x{i}_{j}"))
    rows = []
    for j in range(n_cust):
        rows.append(Row(tuple((xi(i, j), 1.0) for i in range(n_fac)), "==", 1.0, name=f"dem{j}"))
    if form == "weak":
        for i in range(n_fac):
            coeffs = tuple([(xi(i, j), 1.0) for j in range(n_cust)] + [(i, -float(n_cust))])
            rows.append(Row(coeffs, "<=", 0.0, name=f"link{i}"))
    else:
        for i in range(n_fac):
            for j in range(n_cust):
                rows.append(Row(((xi(i, j), 1.0), (i, -1.0)), "<=", 0.0, name=f"link{i}_{j}"))
    obj = [(i, float(f[i])) for i in range(n_fac)]
    obj += [(xi(i, j), float(int(100 * dist[i][j]))) for i in range(n_fac) for j in range(n_cust)]
    planted = [0.0] * len(vs)
    planted[0] = 1.0
    for j in range(n_cust):
        planted[xi(0, j)] = 1.0
    return Model(f"uncap_facility_location_{n_fac}x{n_cust}_{form}", vs, rows, obj, "minimize",
                 meta={"family": "uncap_facility_location", "n_fac": n_fac, "n_cust": n_cust,
                       "form": form, "seed": seed, "planted": planted})


def lot_sizing(T: int = 120, capacitated: bool = True, big_m: float = 0.0,
               seed: int = 42) -> Model:
    """Single-item lot sizing with big-M setups.

        min sum p_t x_t + h_t s_t + f_t y_t
        s.t. s_{t-1} + x_t - s_t = d_t      (s_{-1} = 0)
             x_t <= M_t y_t,  M_t = min(C_t, total demand)   (weak big-M)
             y binary, x, s >= 0
    M uses the TOTAL demand (not the remaining demand) so coefficient
    tightening / (l,S) inequalities have something to do. Probes (c-)MIR via
    path aggregation and big-M tightening. Pochet & Wolsey 2006, "Production
    Planning by Mixed Integer Programming", ch. 7; Barany, Van Roy & Wolsey
    1984, "Uncapacitated lot-sizing: the convex hull of solutions", Math.
    Prog. Study 22. Witness: produce each period's demand (x_t = d_t, y_t = 1),
    feasible because C_t >= d_t.
    big_m > 0 (wide-coefficient variant):
    the setup rows read x_t <= big_m y_t and M_t moves to the column bound
    x_t <= M_t. The integer optimum is unchanged, the LP is weaker, and the
    coefficient range grows to big_m (MIPLIB 2017: 21 of 47 benchmark
    instances have dynamism > 1e3). Bound strengthening recovers the tight
    form. Requires big_m >= max M_t.
    """
    rng = random.Random(seed)
    d = [rng.randint(0, 40) for _ in range(T)]
    D = max(1, sum(d))
    C = [rng.randint(max(d[t], 40), 100) if capacitated else D for t in range(T)]
    M = [float(min(C[t], D)) for t in range(T)]
    if big_m and big_m < max(M):
        raise ValueError(f"big_m {big_m} < max M_t {max(M)}")
    xs = lambda t: t  # noqa: E731
    ss = lambda t: T + t  # noqa: E731
    ys = lambda t: 2 * T + t  # noqa: E731
    vs = ([Var(0.0, M[t] if big_m else INF, name=f"x{t}") for t in range(T)]
          + [Var(0.0, INF, name=f"s{t}") for t in range(T)]
          + _bins(T, "y"))
    rows = []
    for t in range(T):
        coeffs = [(xs(t), 1.0), (ss(t), -1.0)]
        if t > 0:
            coeffs.append((ss(t - 1), 1.0))
        rows.append(Row(tuple(coeffs), "==", float(d[t]), name=f"bal{t}"))
    for t in range(T):
        rows.append(Row(((xs(t), 1.0), (ys(t), -(big_m or M[t]))), "<=", 0.0, name=f"setup{t}"))
    obj = ([(xs(t), float(rng.randint(1, 5))) for t in range(T)]
           + [(ss(t), float(rng.randint(1, 3))) for t in range(T)]
           + [(ys(t), float(rng.randint(30, 150))) for t in range(T)])
    planted = [float(v) for v in d] + [0.0] * T + [1.0] * T
    return Model(f"lot_sizing_{T}{'_bigM' if big_m else ''}", vs, rows, obj, "minimize",
                 meta={"family": "lot_sizing", "T": T, "capacitated": capacitated,
                       "big_m": big_m, "seed": seed, "planted": planted})


def gap(n_agents: int = 10, n_jobs: int = 60, seed: int = 42, infeasible: bool = False) -> Model:
    """Generalized assignment: min c'x s.t. sum_i x_ij = 1, sum_j a_ij x_ij <= b_i.

    Probes knapsack covers on capacity rows + set-partition rows, and diving
    heuristics. Martello & Toth 1981, "An algorithm for the generalized
    assignment problem"; type-D-like data after Chu & Beasley 1997, "A genetic
    algorithm for the generalised assignment problem", C&OR 24: a_ij ~
    U[1,100], c_ij = 111 - a_ij + U[-10,10]. Feasible mode plants a random
    assignment and sets b_i = max(its load, 0.8 * sum_j a_ij / n_agents).
    infeasible=True: n_jobs := 2*n_agents + 1, every a_ij in [34,40], b_i = 100
    -> each agent fits at most 2 jobs (3*34 > 100) -> pigeonhole infeasible,
    while the LP holds (2m+1)*40 <= 100m for m >= 2.
    """
    rng = random.Random(seed)
    m = n_agents
    if infeasible:
        if m < 2:
            raise ValueError("infeasible gap needs n_agents >= 2")
        n_jobs = 2 * m + 1
        a = [[rng.randint(34, 40) for _ in range(n_jobs)] for _ in range(m)]
        b = [100] * m
        planted = None
    else:
        a = [[rng.randint(1, 100) for _ in range(n_jobs)] for _ in range(m)]
        assign = [rng.randrange(m) for _ in range(n_jobs)]
        load = [sum(a[i][j] for j in range(n_jobs) if assign[j] == i) for i in range(m)]
        b = [max(load[i], int(0.8 * sum(a[i]) / m)) for i in range(m)]
        planted = [1 if assign[j] == i else 0 for i in range(m) for j in range(n_jobs)]
    xi = lambda i, j: i * n_jobs + j  # noqa: E731
    vs = [Var(0, 1, VTYPE_BIN, name=f"x{i}_{j}") for i in range(m) for j in range(n_jobs)]
    rows = [Row(tuple((xi(i, j), 1.0) for i in range(m)), "==", 1.0, name=f"assign{j}")
            for j in range(n_jobs)]
    rows += [Row(tuple((xi(i, j), float(a[i][j])) for j in range(n_jobs)), "<=", float(b[i]),
                 name=f"cap{i}") for i in range(m)]
    obj = [(xi(i, j), float(111 - a[i][j] + rng.randint(-10, 10)))
           for i in range(m) for j in range(n_jobs)]
    meta = {"family": "gap", "n_agents": m, "n_jobs": n_jobs, "seed": seed,
            "infeasible": infeasible}
    if planted is not None:
        meta["planted"] = planted
    return Model(f"gap_{m}x{n_jobs}{'_infeas' if infeasible else ''}", vs, rows, obj,
                 "minimize", meta=meta)


def fixed_charge_network(nodes: int = 60, arcs: int = 300, big_m: float = 0.0,
                         var_cost: bool = True, seed: int = 42) -> Model:
    """Single-commodity fixed-charge network flow.

        min sum c_a f_a + F_a y_a
        s.t. out(i) - in(i) = b_i,  0 <= f_a <= u_a y_a,  y binary
    Probes path c-MIR / flow covers and network structure under fixed
    charges. Ortega & Wolsey 2003, "A branch-and-cut algorithm for the
    single-commodity, uncapacitated, fixed-charge network flow problem",
    Networks 41. Feasibility by construction: positive integer flows planted
    on a randomly oriented spanning tree define the balances b; y = 1 exactly
    on the planted-flow arcs.
    big_m > 0: wide-coefficient variant, as in lot_sizing: f_a <= big_m y_a
    with u_a moved to the column bound (same integer optimum, weaker LP).
    var_cost=False: pure fixed-charge objective (c_a = 0, integral F_a only on
    the binaries), the "fixed charge problem" of Hirsch & Dantzig 1968, Naval
    Res. Logist. Q. 15. This is the objective shape of the MIPLIB 2017
    fixed-charge flow instances (beasleyC3, mc11, p200x1188c, sp150x300d:
    obj only on binaries, integral, density 1/2). Same rng draws, so the network is the one var_cost=True builds.
    """
    rng = random.Random(seed)
    arc_list: List[Tuple[int, int]] = []
    arc_set = set()
    planted_f: Dict[Tuple[int, int], int] = {}
    cap: Dict[Tuple[int, int], int] = {}
    for i in range(1, nodes):
        j = rng.randrange(i)
        a = (i, j) if rng.random() < 0.5 else (j, i)
        arc_list.append(a)
        arc_set.add(a)
        cap[a] = rng.randint(10, 40)
        planted_f[a] = rng.randint(1, cap[a])
    guard = 0
    while len(arc_list) < arcs and guard < 20 * arcs:
        guard += 1
        u, v = rng.randrange(nodes), rng.randrange(nodes)
        if u != v and (u, v) not in arc_set:
            arc_set.add((u, v))
            arc_list.append((u, v))
            cap[(u, v)] = rng.randint(10, 40)
    bal = [0] * nodes
    for (u, v), fl in planted_f.items():
        bal[u] += fl
        bal[v] -= fl
    na = len(arc_list)
    if big_m and big_m < max(cap.values()):
        raise ValueError(f"big_m {big_m} < max capacity {max(cap.values())}")
    vs = ([Var(0.0, float(cap[(u, v)]) if big_m else INF, name=f"f{u}_{v}")
           for u, v in arc_list]
          + [Var(0, 1, VTYPE_BIN, name=f"y{u}_{v}") for u, v in arc_list])
    rows = []
    for i in range(nodes):
        coeffs = []
        for k, (u, v) in enumerate(arc_list):
            if u == i:
                coeffs.append((k, 1.0))
            elif v == i:
                coeffs.append((k, -1.0))
        rows.append(Row(tuple(coeffs), "==", float(bal[i]), name=f"bal{i}"))
    for k, a in enumerate(arc_list):
        rows.append(Row(((k, 1.0), (na + k, -float(big_m or cap[a]))), "<=", 0.0,
                        name=f"vub{k}"))
    obj = ([(k, float(rng.randint(1, 10))) for k in range(na)]
           + [(na + k, float(rng.randint(20, 100))) for k in range(na)])
    if not var_cost:
        obj = obj[na:]
    planted = [float(planted_f.get(a, 0)) for a in arc_list]
    planted += [1.0 if a in planted_f else 0.0 for a in arc_list]
    tag = ("_bigM" if big_m else "") + ("" if var_cost else "_fixed")
    return Model(f"fixed_charge_network_{nodes}n{na}a{tag}", vs, rows, obj,
                 "minimize",
                 meta={"family": "fixed_charge_network", "nodes": nodes, "arcs": na,
                       "big_m": big_m, "var_cost": var_cost, "seed": seed,
                       "planted": planted})


def bin_packing(n_items: int = 30, capacity: int = 100, extra_bins: int = 1,
                seed: int = 42) -> Model:
    """Bin packing, assignment (Kantorovich) formulation — highly symmetric.

        min sum_k y_k
        s.t. sum_k x_ik = 1                      (each item packed)
             sum_i w_i x_ik <= C y_k             (bin capacity)
    with K = (first-fit-decreasing bins) + extra_bins, all bins identical ->
    K! symmetric copies of every solution. Probes orbital fixing / orbitopes.
    Martello & Toth 1990, "Knapsack Problems", ch. 8; Kantorovich 1960;
    symmetry: Margot 2010, "Symmetry in integer linear programming".
    Witness: the FFD packing.
    """
    rng = random.Random(seed)
    w = [rng.randint(capacity // 10, capacity // 2 + capacity // 5) for _ in range(n_items)]
    bins: List[List[int]] = []
    loads: List[int] = []
    for i in sorted(range(n_items), key=lambda i: -w[i]):
        for k in range(len(bins)):
            if loads[k] + w[i] <= capacity:
                bins[k].append(i)
                loads[k] += w[i]
                break
        else:
            bins.append([i])
            loads.append(w[i])
    K = len(bins) + extra_bins
    xi = lambda i, k: i * K + k  # noqa: E731
    vs = [Var(0, 1, VTYPE_BIN, name=f"x{i}_{k}") for i in range(n_items) for k in range(K)]
    yi = lambda k: n_items * K + k  # noqa: E731
    vs += _bins(K, "y")
    rows = [Row(tuple((xi(i, k), 1.0) for k in range(K)), "==", 1.0, name=f"item{i}")
            for i in range(n_items)]
    rows += [Row(tuple([(xi(i, k), float(w[i])) for i in range(n_items)]
                       + [(yi(k), -float(capacity))]), "<=", 0.0, name=f"bin{k}")
             for k in range(K)]
    obj = [(yi(k), 1.0) for k in range(K)]
    planted = [0.0] * len(vs)
    for k, items in enumerate(bins):
        planted[yi(k)] = 1.0
        for i in items:
            planted[xi(i, k)] = 1.0
    return Model(f"bin_packing_{n_items}", vs, rows, obj, "minimize",
                 meta={"family": "bin_packing", "n_items": n_items, "capacity": capacity,
                       "bins": K, "seed": seed, "planted": planted,
                       "planted_obj": float(len(bins))})


def _random_graph(rng, n, density):
    return [(u, v) for u in range(n) for v in range(u + 1, n) if rng.random() < density]


def graph_coloring(n_nodes: int = 20, density: float = 0.3, max_edges: int = 0,
                   seed: int = 42) -> Model:
    """Vertex coloring, assignment formulation — symmetric in the colors.

        min sum_c w_c
        s.t. sum_c x_vc = 1                  (each vertex colored)
             x_uc + x_vc <= w_c              (each edge, each color)
    K = greedy-coloring colors, so a witness always exists; max_edges > 0
    keeps a random subset of that many edges (bounds the row count for the
    `small` size). Probes symmetry
    handling (color permutations) and clique cuts. Mendez-Diaz & Zabala 2006,
    "A branch-and-cut algorithm for graph coloring", Discrete Appl. Math. 154.
    """
    rng = random.Random(seed)
    edges = _random_graph(rng, n_nodes, density)
    if 0 < max_edges < len(edges):
        edges = sorted(rng.sample(edges, max_edges))
    adj = [set() for _ in range(n_nodes)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    color = [-1] * n_nodes
    for v in sorted(range(n_nodes), key=lambda v: -len(adj[v])):
        used = {color[u] for u in adj[v]}
        c = 0
        while c in used:
            c += 1
        color[v] = c
    K = max(color) + 1
    xi = lambda v, c: v * K + c  # noqa: E731
    wi = lambda c: n_nodes * K + c  # noqa: E731
    vs = [Var(0, 1, VTYPE_BIN, name=f"x{v}_{c}") for v in range(n_nodes) for c in range(K)]
    vs += _bins(K, "w")
    rows = [Row(tuple((xi(v, c), 1.0) for c in range(K)), "==", 1.0, name=f"col{v}")
            for v in range(n_nodes)]
    for u, v in edges:
        for c in range(K):
            rows.append(Row(((xi(u, c), 1.0), (xi(v, c), 1.0), (wi(c), -1.0)), "<=", 0.0,
                            name=f"e{u}_{v}_{c}"))
    obj = [(wi(c), 1.0) for c in range(K)]
    planted = [0.0] * len(vs)
    for v in range(n_nodes):
        planted[xi(v, color[v])] = 1.0
    for c in range(K):
        planted[wi(c)] = 1.0
    return Model(f"graph_coloring_{n_nodes}", vs, rows, obj, "minimize",
                 meta={"family": "graph_coloring", "n_nodes": n_nodes, "density": density,
                       "edges": len(edges), "colors": K, "seed": seed, "planted": planted,
                       "planted_obj": float(K)})


def max_indep_set(n_nodes: int = 150, density: float = 0.05, weighted: bool = False,
                  seed: int = 42) -> Model:
    """Maximum (weighted) independent set, edge formulation x_u + x_v <= 1.

    The edge formulation's LP is half-integral (x = 1/2 everywhere) -> clique
    extraction from the edge rows + clique cuts carry the bound. Nemhauser &
    Trotter 1974, "Properties of vertex packing and independence system
    polyhedra", Math. Prog. 6; Padberg 1973 (clique facets). Witness: greedy
    independent set.
    """
    rng = random.Random(seed)
    edges = _random_graph(rng, n_nodes, density)
    rows = [Row(((u, 1.0), (v, 1.0)), "<=", 1.0, name=f"e{u}_{v}") for u, v in edges]
    obj = [(v, float(rng.randint(1, 10)) if weighted else 1.0) for v in range(n_nodes)]
    adj = [set() for _ in range(n_nodes)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    chosen: set = set()
    for v in sorted(range(n_nodes), key=lambda v: len(adj[v])):
        if not adj[v] & chosen:
            chosen.add(v)
    planted = [1 if v in chosen else 0 for v in range(n_nodes)]
    return Model(f"max_indep_set_{n_nodes}", _bins(n_nodes), rows, obj, "maximize",
                 meta={"family": "max_indep_set", "n_nodes": n_nodes, "density": density,
                       "edges": len(edges), "weighted": weighted, "seed": seed,
                       "planted": planted})


def market_split(m: int = 4, n: int = 0, seed: int = 42) -> Model:
    """Market split (Cornuejols & Dawande 1999), optimization version.

        min sum_i (sp_i + sm_i)
        s.t. sum_j a_ij x_j + sp_i - sm_i = d_i,  x binary, sp, sm >= 0
    a_ij ~ U[0,99], d_i = floor(sum_j a_ij / 2), n = 10(m-1) by default.
    Tiny but pathological for LP-based B&B (cuts barely help; branching
    only). Cornuejols & Dawande 1999, "A class of hard small 0-1 programs",
    INFORMS J. Computing 11. Witness: x = 0, sp = d (objective sum d).
    """
    rng = random.Random(seed)
    if n <= 0:
        n = 10 * (m - 1)
    a = [[rng.randint(0, 99) for _ in range(n)] for _ in range(m)]
    d = [sum(a[i]) // 2 for i in range(m)]
    vs = _bins(n) + [Var(0.0, INF, name=f"sp{i}") for i in range(m)] \
        + [Var(0.0, INF, name=f"sm{i}") for i in range(m)]
    rows = [Row(tuple([(j, float(a[i][j])) for j in range(n) if a[i][j]]
                      + [(n + i, 1.0), (n + m + i, -1.0)]), "==", float(d[i]), name=f"split{i}")
            for i in range(m)]
    obj = [(n + i, 1.0) for i in range(m)] + [(n + m + i, 1.0) for i in range(m)]
    planted = [0.0] * n + [float(v) for v in d] + [0.0] * m
    return Model(f"market_split_{m}x{n}", vs, rows, obj, "minimize",
                 meta={"family": "market_split", "m": m, "n": n, "seed": seed,
                       "planted": planted, "planted_obj": float(sum(d))})


def jobshop_bigM(n_jobs: int = 10, n_machines: int = 5, big_m: float = 0.0, p_max: int = 10,
                 seed: int = 42) -> Model:
    """Job-shop makespan, Manne's disjunctive big-M formulation.

        min Cmax
        s.t. S_{j,o+1} >= S_{j,o} + p_{j,o}                    (job routing)
             Cmax >= S_{j,last} + p_{j,last}
             S_{i,k} >= S_{j,k} + p_{j,k} - M z_{ijk}          (i after j, or)
             S_{j,k} >= S_{i,k} + p_{i,k} - M (1 - z_{ijk})    (j after i)
    with M = sum of all processing times. The LP relaxation is near-useless
    (z = 1/2 kills every disjunction) -> branching- and propagation-
    dominated. Manne 1960, "On the job-shop scheduling problem", Oper. Res. 8.
    Witness: run the jobs one after another (no overlap anywhere).
    Wide-coefficient variant: p_max
    widens the processing times p ~ U[1, p_max], and big_m > 0 replaces the
    horizon M = sum p by a modeller's "large constant" (must be >= sum p).
    That practice is what drives MIPLIB coefficient ranges to 1e6 and beyond
    (glass4: 8.4e6).
    """
    rng = random.Random(seed)
    routes = []
    for _ in range(n_jobs):
        r = list(range(n_machines))
        rng.shuffle(r)
        routes.append(r)
    p = [[rng.randint(1, p_max) for _ in range(n_machines)] for _ in range(n_jobs)]  # p[j][k]
    M = float(sum(map(sum, p)))
    if big_m:
        if big_m < M:
            raise ValueError(f"big_m {big_m} < horizon {M}")
        M = float(big_m)
    si = lambda j, k: j * n_machines + k  # noqa: E731  start of job j on machine k
    cmax = n_jobs * n_machines
    vs = [Var(0.0, INF, name=f"S{j}_{k}") for j in range(n_jobs) for k in range(n_machines)]
    vs.append(Var(0.0, INF, name="Cmax"))
    pairs = [(i, j, k) for k in range(n_machines) for i in range(n_jobs)
             for j in range(i + 1, n_jobs)]
    zi = {pk: cmax + 1 + t for t, pk in enumerate(pairs)}
    vs += [Var(0, 1, VTYPE_BIN, name=f"z{i}_{j}_{k}") for i, j, k in pairs]
    rows = []
    for j in range(n_jobs):
        r = routes[j]
        for o in range(n_machines - 1):
            rows.append(Row(((si(j, r[o + 1]), 1.0), (si(j, r[o]), -1.0)), ">=",
                            float(p[j][r[o]]), name=f"prec{j}_{o}"))
        rows.append(Row(((cmax, 1.0), (si(j, r[-1]), -1.0)), ">=", float(p[j][r[-1]]),
                        name=f"mk{j}"))
    for (i, j, k) in pairs:
        z = zi[(i, j, k)]
        # S_ik - S_jk + M z >= p_jk ; S_jk - S_ik - M z >= p_ik - M
        rows.append(Row(((si(i, k), 1.0), (si(j, k), -1.0), (z, M)), ">=", float(p[j][k]),
                        name=f"dA{i}_{j}_{k}"))
        rows.append(Row(((si(j, k), 1.0), (si(i, k), -1.0), (z, -M)), ">=", float(p[i][k]) - M,
                        name=f"dB{i}_{j}_{k}"))
    obj = [(cmax, 1.0)]
    # z_ijk = 1 <=> job j follows job i on machine k (dB active, dA relaxed
    # by M). Serial witness: job 0 fully, then job 1, ... so every i < j has
    # j after i -> z = 1 everywhere; dA then reads S_jk + p_jk - S_ik <= M,
    # true because every job finishes by sum(p) = M.
    planted = [0.0] * (cmax + 1) + [1.0] * len(pairs)
    t = 0.0
    for j in range(n_jobs):
        for k in routes[j]:
            planted[si(j, k)] = t
            t += p[j][k]
    planted[cmax] = t
    tag = f"_M{big_m:g}" if big_m else ""
    return Model(f"jobshop_bigM_{n_jobs}x{n_machines}{tag}", vs, rows, obj, "minimize",
                 meta={"family": "jobshop_bigM", "n_jobs": n_jobs, "n_machines": n_machines,
                       "seed": seed, "big_m": M, "planted": planted})


# ---------------------------------------------------------------------------
# Mixed-structure families. The families above are each one or two MIPLIB 2017 row
# types. Real instances mix 1-5 types (median 3), have wide coefficient ranges,
# general integers and several components. Each family below combines several
# of the 17 row classes of Gleixner et al., "MIPLIB 2017", Math. Prog. Comp. 13
# (2021) §4.4; the docstrings list the classes each row group lands in, as
# `rusolve-cli classify` reports them (default mode).


def _shortest_path(adj: Dict[int, List[Tuple[int, int]]], s: int, t: int) -> List[int]:
    """BFS path s -> t as a list of arc indices (adj[u] = [(v, arc), ...])."""
    prev: Dict[int, Tuple[int, int]] = {s: (-1, -1)}
    frontier = [s]
    while frontier and t not in prev:
        nxt = []
        for u in frontier:
            for v, a in adj.get(u, ()):
                if v not in prev:
                    prev[v] = (u, a)
                    nxt.append(v)
        frontier = nxt
    path, v = [], t
    while v != s:
        u, a = prev[v]
        path.append(a)
        v = u
    return path[::-1]


def network_design(nodes: int = 12, arcs: int = 36, commodities: int = 6, strong: bool = True,
                   cost_scale: float = 1.0, var_cost: bool = True, modules: int = 0,
                   seed: int = 42) -> Model:
    """Capacitated multicommodity fixed-charge network design (MCND) with side rows.

        min  sum_a c_a sum_k x_a^k + sum_a f_a y_a
        s.t. out_k(i) - in_k(i) = d_k [i = o_k] - d_k [i = t_k]  (flow, mixed binary)
             sum_k x_a^k <= u_a y_a                (arc capacity, mixed binary)
             x_a^k <= min(d_k, u_a) y_a            (strong linking, variable bound)
             sum_a f_a y_a <= B                    (design budget, knapsack)
             sum_{a in out(i)} y_a <= deg_i        (degree limit, invariant knapsack)
             x >= 0, y binary
    Formulation: Gendron, Crainic & Frangioni 1999, "Multicommodity capacitated
    network design", in Telecommunications Network Planning (Kluwer); the
    strong linking rows are what separates their strong from their weak
    relaxation. Budget and degree rows are the usual side constraints in
    survivable/telecom design (e.g. Magnanti & Wong 1984, "Network design and
    transportation planning: models and algorithms", Transp. Sci. 18).
    strong=False drops the A*K linking rows (weak form).
    cost_scale multiplies the fixed costs f_a, widening the objective range.
    var_cost=False drops the flow costs c_a (objective on the design variables
    only, as in the fixed-charge variant of fixed_charge_network).
    modules=K > 0 is network loading: y_a becomes an integer number of
    capacity modules in [0, K], each of capacity u_a, so the capacity row reads
    sum_k x_a^k <= u_a y_a with y general integer (Magnanti, Mirchandani &
    Vachani 1995, "Modeling and solving the two-facility capacitated network
    loading problem", Oper. Res. 43; Atamtürk & Rajan 2002, "On splittable and
    unsplittable flow capacitated network design arc-set polyhedra", Math.
    Prog. 92). Budget and degree rows then become integer knapsacks. The
    planted design installs ceil(load / u_a) modules.
    Feasibility by construction: a Hamiltonian cycle makes the graph strongly
    connected. Each commodity is routed on a BFS path, capacities are raised to
    the planted load, B = 1.1 x the planted design cost, and deg_i = planted
    out-degree + 1 (a row is kept only when it cuts something).
    """
    rng = random.Random(seed)
    nodes = max(nodes, 3)
    order = list(range(nodes))
    rng.shuffle(order)
    arc_list: List[Tuple[int, int]] = [(order[i], order[(i + 1) % nodes]) for i in range(nodes)]
    arc_set = set(arc_list)
    guard = 0
    while len(arc_list) < arcs and guard < 50 * arcs:
        guard += 1
        u, v = rng.randrange(nodes), rng.randrange(nodes)
        if u != v and (u, v) not in arc_set:
            arc_set.add((u, v))
            arc_list.append((u, v))
    na, K = len(arc_list), commodities
    adj: Dict[int, List[Tuple[int, int]]] = {}
    for a in sorted(range(na), key=lambda _: rng.random()):
        u, v = arc_list[a]
        adj.setdefault(u, []).append((v, a))
    od, dem, paths = [], [], []
    for _ in range(K):
        o, t = rng.sample(range(nodes), 2)
        od.append((o, t))
        dem.append(rng.randint(5, 25))
        paths.append(_shortest_path(adj, o, t))
    load = [0] * na
    for k in range(K):
        for a in paths[k]:
            load[a] += dem[k]
    cap = [max(rng.randint(15, 60), load[a]) for a in range(na)]
    c = [rng.randint(1, 10) for _ in range(na)]
    f = [int(round(rng.randint(40, 200) * cost_scale)) for _ in range(na)]
    if modules:
        # module capacity: the planted load fits in `modules` modules
        cap = [max(rng.randint(8, 30), -(-load[a] // modules)) for a in range(na)]
    ny = [(-(-load[a] // cap[a]) if load[a] else 0) for a in range(na)]
    xi = lambda a, k: a * K + k  # noqa: E731
    yi = lambda a: na * K + a  # noqa: E731
    vs = [Var(0.0, INF, name=f"x{a}_{k}") for a in range(na) for k in range(K)]
    if modules:
        vs += [Var(0, float(modules), VTYPE_INT, name=f"y{a}") for a in range(na)]
    else:
        vs += [Var(0, 1, VTYPE_BIN, name=f"y{a}") for a in range(na)]
    rows: List[Row] = []
    for k in range(K):
        o, t = od[k]
        for i in range(nodes):
            coeffs = [(xi(a, k), 1.0 if u == i else -1.0)
                      for a, (u, v) in enumerate(arc_list) if i in (u, v)]
            rhs = float(dem[k] if i == o else (-dem[k] if i == t else 0))
            rows.append(Row(tuple(coeffs), "==", rhs, name=f"flow{k}_{i}"))
    for a in range(na):
        rows.append(Row(tuple([(xi(a, k), 1.0) for k in range(K)] + [(yi(a), -float(cap[a]))]),
                        "<=", 0.0, name=f"cap{a}"))
    if strong:
        for a in range(na):
            for k in range(K):
                lim = min(dem[k], cap[a] * max(1, modules))
                rows.append(Row(((xi(a, k), 1.0), (yi(a), -float(lim))), "<=",
                                0.0, name=f"link{a}_{k}"))
    used = {a for p in paths for a in p}
    budget = int(1.1 * sum(f[a] * ny[a] for a in used)) + 1
    rows.append(Row(tuple((yi(a), float(f[a])) for a in range(na)), "<=", float(budget),
                    name="budget"))
    for i in range(nodes):
        out = [a for a, (u, _) in enumerate(arc_list) if u == i]
        deg = sum(ny[a] for a in out if a in used) + 1
        if len(out) * max(1, modules) > deg:
            rows.append(Row(tuple((yi(a), 1.0) for a in out), "<=", float(deg), name=f"deg{i}"))
    obj = ([(xi(a, k), float(c[a])) for a in range(na) for k in range(K)] if var_cost else [])
    obj += [(yi(a), float(f[a])) for a in range(na)]
    planted = [0.0] * len(vs)
    for k in range(K):
        for a in paths[k]:
            planted[xi(a, k)] += dem[k]
    for a in used:
        planted[yi(a)] = float(ny[a])
    tag = ("" if strong else "_weak") + ("" if var_cost else "_fixed") + \
        (f"_mod{modules}" if modules else "")
    return Model(f"network_design_{nodes}n{na}a{K}k{tag}", vs, rows, obj,
                 "minimize",
                 meta={"family": "network_design", "nodes": nodes, "arcs": na, "commodities": K,
                       "strong": strong, "cost_scale": cost_scale, "var_cost": var_cost,
                       "modules": modules, "seed": seed, "planted": planted})


def facility_mixed(n_fac: int = 8, n_cust: int = 20, hubs: int = 2, p_open: int = 0,
                   module_cap: int = 40, max_modules: int = 3, link_frac: float = 0.3,
                   cover_frac: float = 0.3, seed: int = 42) -> Model:
    """Hierarchical modular capacitated facility location with side rows.

    y_i binary (open), z_i integer (capacity modules), x_ij in [0, 1] (share of
    customer j's demand d_j served by i):
        min  sum f_i y_i + sum g_i z_i + sum c_ij d_j x_ij
        s.t. sum_i x_ij = 1                         (demand, mixed binary)
             sum_j d_j x_ij - Q z_i <= 0            (modular capacity, general linear)
             z_i - K_i y_i <= 0                     (modules need an open site, variable bound)
             x_ij - y_i <= 0, nearest pairs only    (strong linking, variable bound)
             y_s - y_h(s) <= 0                      (satellite needs its hub, precedence)
             sum_i y_i = p                          (p sites, cardinality)
             sum f_i y_i + sum g_i z_i <= B         (budget, integer knapsack)
             sum_{i near j} y_i >= 1, some j        (service coverage, set covering)
    Sources: modular capacities after Correia & Captivo 2003, "A Lagrangean
    heuristic for a modular capacitated location problem", Ann. Oper. Res. 122;
    hub/satellite precedence after Sahin & Sural 2007, "A review of
    hierarchical facility location models", Comput. Oper. Res. 34; the
    p-facility cardinality row as in the p-median (ReVelle & Swain 1970,
    "Central facilities location", Geogr. Anal. 2); coverage rows as in
    Church & ReVelle 1974, "The maximal covering location problem", Pap. Reg.
    Sci. Assoc. 32.
    Witness: open all hubs plus random satellites (p in total), assign each
    customer to its nearest open site, z_i = ceil(load_i / Q). B = 1.1 x the
    planted fixed cost. A covering row for customer j lists the sites within
    1.3 x its planted distance (which includes the planted site).
    """
    rng = random.Random(seed)
    hubs = max(1, min(hubs, n_fac - 1))
    dist = _facility_data(rng, n_fac, n_cust)
    d = [rng.randint(5, 30) for _ in range(n_cust)]
    parent = [-1] * hubs + [rng.randrange(hubs) for _ in range(n_fac - hubs)]
    p = p_open if p_open > 0 else max(hubs + 1, n_fac // 2)
    p = min(max(p, hubs), n_fac)
    open_ = set(range(hubs)) | set(rng.sample(range(hubs, n_fac), p - hubs))
    assign = [min(open_, key=lambda i: dist[i][j]) for j in range(n_cust)]
    load = [sum(d[j] for j in range(n_cust) if assign[j] == i) for i in range(n_fac)]
    Q = module_cap
    zp = [-(-load[i] // Q) for i in range(n_fac)]
    Kmax = [max(max_modules, zp[i]) for i in range(n_fac)]
    f = [rng.randint(300, 900) for _ in range(n_fac)]
    g = [rng.randint(60, 200) for _ in range(n_fac)]
    yi = lambda i: i  # noqa: E731
    zi = lambda i: n_fac + i  # noqa: E731
    xi = lambda i, j: 2 * n_fac + i * n_cust + j  # noqa: E731
    vs = _bins(n_fac, "y") + [Var(0, float(Kmax[i]), VTYPE_INT, name=f"z{i}")
                              for i in range(n_fac)]
    vs += [Var(0.0, 1.0, name=f"x{i}_{j}") for i in range(n_fac) for j in range(n_cust)]
    rows: List[Row] = []
    for j in range(n_cust):
        rows.append(Row(tuple((xi(i, j), 1.0) for i in range(n_fac)), "==", 1.0, name=f"dem{j}"))
    for i in range(n_fac):
        rows.append(Row(tuple([(xi(i, j), float(d[j])) for j in range(n_cust)]
                              + [(zi(i), -float(Q))]), "<=", 0.0, name=f"cap{i}"))
    for i in range(n_fac):
        rows.append(Row(((zi(i), 1.0), (yi(i), -float(Kmax[i]))), "<=", 0.0, name=f"mod{i}"))
    pairs = sorted(((dist[i][j], i, j) for i in range(n_fac) for j in range(n_cust)))
    for _, i, j in pairs[:int(link_frac * len(pairs))]:
        rows.append(Row(((xi(i, j), 1.0), (yi(i), -1.0)), "<=", 0.0, name=f"link{i}_{j}"))
    for s in range(hubs, n_fac):
        rows.append(Row(((yi(s), 1.0), (yi(parent[s]), -1.0)), "<=", 0.0, name=f"hier{s}"))
    rows.append(Row(tuple((yi(i), 1.0) for i in range(n_fac)), "==", float(p), name="card"))
    budget = int(1.1 * sum(f[i] + g[i] * zp[i] for i in open_)) + 1
    rows.append(Row(tuple([(yi(i), float(f[i])) for i in range(n_fac)]
                          + [(zi(i), float(g[i])) for i in range(n_fac)]), "<=", float(budget),
                    name="budget"))
    for j in sorted(rng.sample(range(n_cust), int(cover_frac * n_cust))):
        near = [i for i in range(n_fac) if dist[i][j] <= 1.3 * dist[assign[j]][j] + 1e-12]
        if len(near) >= 2:
            rows.append(Row(tuple((yi(i), 1.0) for i in near), ">=", 1.0, name=f"cov{j}"))
    obj = [(yi(i), float(f[i])) for i in range(n_fac)] + [(zi(i), float(g[i]))
                                                         for i in range(n_fac)]
    obj += [(xi(i, j), float(int(100 * dist[i][j] * d[j]) / 10))
            for i in range(n_fac) for j in range(n_cust)]
    planted = [0.0] * len(vs)
    for i in open_:
        planted[yi(i)] = 1.0
        planted[zi(i)] = float(zp[i])
    for j in range(n_cust):
        planted[xi(assign[j], j)] = 1.0
    return Model(f"facility_mixed_{n_fac}x{n_cust}", vs, rows, obj, "minimize",
                 meta={"family": "facility_mixed", "n_fac": n_fac, "n_cust": n_cust,
                       "hubs": hubs, "p": p, "seed": seed, "planted": planted})


def budget_max_cover(n_sets: int = 60, n_elems: int = 120, density: float = 0.05, k: int = 8,
                     budget_frac: float = 0.6, conflicts: int = 10, seed: int = 42) -> Model:
    """Budgeted maximum coverage with a cardinality limit and conflicts.

        max  sum_i w_i z_i
        s.t. z_i - sum_{j covers i} x_j <= 0   (coverage, invariant knapsack after
                                                complementing; set packing if one set)
             sum_j x_j <= k                    (cardinality budget, invariant knapsack)
             sum_j c_j x_j <= B                (cost budget, knapsack)
             x_a + x_b <= 1                    (conflicting pairs, variable bound)
             x, z binary
    Khuller, Moss & Naor 1999, "The budgeted maximum coverage problem", Inf.
    Process. Lett. 70; maximum coverage under a cardinality constraint:
    Nemhauser, Wolsey & Fisher 1978, "An analysis of approximations for
    maximizing submodular set functions", Math. Prog. 14. This family puts the
    cardinality/budget rows of the corpus recommendations next to covering-type
    structure. B = budget_frac x the cost of k average sets. Witness: greedy
    over a random set order within k and B. Conflict pairs are never both in
    the witness.
    """
    rng = random.Random(seed)
    covers: List[List[int]] = []
    for _ in range(n_sets):
        s = [i for i in range(n_elems) if rng.random() < density] or [rng.randrange(n_elems)]
        covers.append(s)
    by_elem: List[List[int]] = [[] for _ in range(n_elems)]
    for j, s in enumerate(covers):
        for i in s:
            by_elem[i].append(j)
    cost = [rng.randint(10, 100) for _ in range(n_sets)]
    w = [rng.randint(1, 50) for _ in range(n_elems)]
    B = int(budget_frac * k * sum(cost) / n_sets)
    chosen: List[int] = []
    spent = 0
    for j in sorted(range(n_sets), key=lambda _: rng.random()):
        if len(chosen) < k and spent + cost[j] <= B:
            chosen.append(j)
            spent += cost[j]
    zi = lambda i: n_sets + i  # noqa: E731
    vs = _bins(n_sets, "x") + _bins(n_elems, "z")
    rows: List[Row] = []
    for i in range(n_elems):
        if by_elem[i]:
            rows.append(Row(tuple([(zi(i), 1.0)] + [(j, -1.0) for j in by_elem[i]]), "<=", 0.0,
                            name=f"cov{i}"))
    rows.append(Row(tuple((j, 1.0) for j in range(n_sets)), "<=", float(k), name="card"))
    rows.append(Row(tuple((j, float(cost[j])) for j in range(n_sets)), "<=", float(B),
                    name="budget"))
    chosen_set = set(chosen)
    seen = set()
    guard = 0
    while len(seen) < conflicts and guard < 100 * max(1, conflicts):
        guard += 1
        a, b = sorted(rng.sample(range(n_sets), 2))
        if (a, b) not in seen and not (a in chosen_set and b in chosen_set):
            seen.add((a, b))
            rows.append(Row(((a, 1.0), (b, 1.0)), "<=", 1.0, name=f"conf{a}_{b}"))
    obj = [(zi(i), float(w[i])) for i in range(n_elems) if by_elem[i]]
    planted = [0.0] * len(vs)
    for j in chosen:
        planted[j] = 1.0
        for i in covers[j]:
            planted[zi(i)] = 1.0
    return Model(f"budget_max_cover_{n_sets}x{n_elems}", vs, rows, obj, "maximize",
                 meta={"family": "budget_max_cover", "n_sets": n_sets, "n_elems": n_elems,
                       "k": k, "budget": B, "seed": seed, "planted": planted})


def int_production(products: int = 6, periods: int = 8, wide: bool = True,
                   seed: int = 42) -> Model:
    """Multi-item production planning with integer batches and an integer workforce.

    Per product p and period t: n_pt integer batches (size B_p), s_pt >= 0
    inventory, y_pt binary setup. Per period: w_t integer workers, o_t >= 0
    overtime hours.
        min  sum k_p n_pt + h_p s_pt + F_p y_pt + sum W w_t + V o_t
        s.t. s_p,t-1 + B_p n_pt - s_pt = d_pt         (balance, general linear)
             n_pt - U_p y_pt <= 0                     (setup, variable bound; U_p big-M)
             sum_p a_p n_pt - H w_t - o_t <= 0        (labour hours, general linear)
             o_t - 0.25 H w_t <= 0                    (overtime cap, general linear)
             w_t - w_t-1 <= hire,  w_t-1 - w_t <= fire (precedence)
             sum_p y_pt <= S                          (setup slots, invariant knapsack)
    Sources: big-bucket multi-item lot sizing (Pochet & Wolsey 2006,
    "Production Planning by Mixed Integer Programming", Springer, ch. 9-11)
    combined with the workforce and overtime rows of the "factory planning"
    and "manpower planning" models in Williams 2013, "Model Building in
    Mathematical Programming", 5th ed., Wiley.
    wide=True (default) gives the MIPLIB-style numerics this corpus lacked:
    general integers with bounds up to 2e4 (MIPLIB's wide_integer > 1e4
    domain), batch sizes and hours spanning 10^-1.5..10^2.5, and costs over
    four decades. wide=False keeps every coefficient in [1, 100].
    Witness: n_pt = ceil(d_pt / B_p), so the surplus carries into inventory;
    w_t = ceil(hours_t / H), o_t = 0; hire/fire/S are the planted maxima + 1.
    """
    rng = random.Random(seed)
    P, T = products, periods

    def logu(lo, hi):
        return 10 ** rng.uniform(lo, hi)

    if wide:
        Bsz = [round(logu(-0.5, 2.5), 2) for _ in range(P)]
        a = [round(logu(-1.5, 1.5), 3) for _ in range(P)]
        U = [rng.randint(1000, 20000) for _ in range(P)]
        F = [round(logu(2, 4)) for _ in range(P)]
        W = rng.randint(500, 2000)
    else:
        Bsz = [float(rng.randint(5, 20)) for _ in range(P)]
        a = [float(rng.randint(1, 5)) for _ in range(P)]
        U = [50] * P
        F = [rng.randint(20, 100) for _ in range(P)]
        W = rng.randint(20, 60)
    H = 40.0
    hcost = [round(rng.uniform(0.1, 2.0), 2) for _ in range(P)]
    kcost = [round(Bsz[p] * rng.uniform(1, 3), 2) for p in range(P)]
    V = round(1.5 * W / H, 2)
    dem = [[round(rng.uniform(0, 6) * Bsz[p], 2) if rng.random() < 0.8 else 0.0
            for _ in range(T)] for p in range(P)]
    ni = lambda p, t: (p * T + t)  # noqa: E731
    si = lambda p, t: P * T + p * T + t  # noqa: E731
    yi = lambda p, t: 2 * P * T + p * T + t  # noqa: E731
    wi = lambda t: 3 * P * T + t  # noqa: E731
    oi = lambda t: 3 * P * T + T + t  # noqa: E731
    # witness
    n0 = [[-(-dem[p][t] // Bsz[p]) for t in range(T)] for p in range(P)]
    for p in range(P):
        for t in range(T):
            n0[p][t] = float(min(int(n0[p][t]), U[p]))
    s0 = [[0.0] * T for _ in range(P)]
    for p in range(P):
        inv = 0.0
        for t in range(T):
            inv = inv + Bsz[p] * n0[p][t] - dem[p][t]
            s0[p][t] = inv
    hours = [sum(a[p] * n0[p][t] for p in range(P)) for t in range(T)]
    w0 = [float(max(1, -(-hours[t] // H))) for t in range(T)]
    hire = max([w0[t] - w0[t - 1] for t in range(1, T)] + [0.0]) + 1
    fire = max([w0[t - 1] - w0[t] for t in range(1, T)] + [0.0]) + 1
    slots = max(2, max(sum(1 for p in range(P) if n0[p][t] > 0) for t in range(T)) + 1)
    wmax = float(max(w0) * 2 + 5)
    vs = [Var(0, float(U[p]), VTYPE_INT, name=f"n{p}_{t}") for p in range(P) for t in range(T)]
    vs += [Var(0.0, INF, name=f"s{p}_{t}") for p in range(P) for t in range(T)]
    vs += [Var(0, 1, VTYPE_BIN, name=f"y{p}_{t}") for p in range(P) for t in range(T)]
    vs += [Var(0, wmax, VTYPE_INT, name=f"w{t}") for t in range(T)]
    vs += [Var(0.0, INF, name=f"o{t}") for t in range(T)]
    rows: List[Row] = []
    for p in range(P):
        for t in range(T):
            coeffs = [(ni(p, t), Bsz[p]), (si(p, t), -1.0)]
            if t > 0:
                coeffs.append((si(p, t - 1), 1.0))
            rows.append(Row(tuple(coeffs), "==", dem[p][t], name=f"bal{p}_{t}"))
    for p in range(P):
        for t in range(T):
            rows.append(Row(((ni(p, t), 1.0), (yi(p, t), -float(U[p]))), "<=", 0.0,
                            name=f"setup{p}_{t}"))
    for t in range(T):
        rows.append(Row(tuple([(ni(p, t), a[p]) for p in range(P)]
                              + [(wi(t), -H), (oi(t), -1.0)]), "<=", 0.0, name=f"hours{t}"))
        rows.append(Row(((oi(t), 1.0), (wi(t), -0.25 * H)), "<=", 0.0, name=f"ot{t}"))
    for t in range(1, T):
        rows.append(Row(((wi(t), 1.0), (wi(t - 1), -1.0)), "<=", hire, name=f"hire{t}"))
        rows.append(Row(((wi(t - 1), 1.0), (wi(t), -1.0)), "<=", fire, name=f"fire{t}"))
    for t in range(T):
        rows.append(Row(tuple((yi(p, t), 1.0) for p in range(P)), "<=", float(slots),
                        name=f"slots{t}"))
    obj = [(ni(p, t), kcost[p]) for p in range(P) for t in range(T)]
    obj += [(si(p, t), hcost[p]) for p in range(P) for t in range(T)]
    obj += [(yi(p, t), float(F[p])) for p in range(P) for t in range(T)]
    obj += [(wi(t), float(W)) for t in range(T)] + [(oi(t), V) for t in range(T)]
    planted = [0.0] * len(vs)
    for p in range(P):
        for t in range(T):
            planted[ni(p, t)] = n0[p][t]
            planted[si(p, t)] = s0[p][t]
            planted[yi(p, t)] = 1.0 if n0[p][t] > 0 else 0.0
    for t in range(T):
        planted[wi(t)] = w0[t]
    return Model(f"int_production_{P}x{T}{'' if wide else '_narrow'}", vs, rows, obj, "minimize",
                 meta={"family": "int_production", "products": P, "periods": T, "wide": wide,
                       "seed": seed, "planted": planted})


def general_ip(m: int = 20, n: int = 30, maximize: bool = True, seed: int = 42) -> Model:
    """Dense all-general-integer program with fractional coefficients and no upper bounds.

        max c'x  s.t.  A x <= b   (maximize=True)
        min c'x  s.t.  A x >= b   (maximize=False)
        x >= 0 integer, no upper bounds
    A is dense and positive with log-uniform entries in [10^-2.5, 10^1.25]
    (4 decimals), c has 2 decimals in [1, 50]. Every row is general linear in
    the MIPLIB 2017 classification (integer columns, fractional coefficients),
    no column is boxed. This is the shape of the MIPLIB 2017 gen-ip instances
    (gen-ip002, gen-ip054: 24-27 rows x 30-41 integer columns, coefficients in
    [0.003, 18]; Gleixner et al. 2021, "MIPLIB 2017", Math. Prog. Comp. 13,
    as `rusolve-cli classify` reports them). Dense
    rows with many near-equal ratios are what makes small pure IPs hard for
    LP-based branch-and-bound (Cornuéjols & Dawande 1999, "A class of hard
    small 0-1 programs", INFORMS J. Comput. 11, the same effect with binaries).
    A positive dense A bounds every column (<= rows) or the objective
    (c > 0 with >= rows), so the IP has an optimum.
    Witness: x_hat ~ U{0..10}; b = A x_hat + U[0, 5] (<=) or max(0, A x_hat - U[0, 5]) (>=).
    """
    rng = random.Random(seed)
    A = [[round(10 ** rng.uniform(-2.5, 1.25), 4) for _ in range(n)] for _ in range(m)]
    xh = [rng.randint(0, 10) for _ in range(n)]
    act = [sum(A[i][j] * xh[j] for j in range(n)) for i in range(m)]
    vs = [Var(0.0, INF, VTYPE_INT, name=f"x{j}") for j in range(n)]
    rows = []
    for i in range(m):
        if maximize:
            rows.append(Row(tuple((j, A[i][j]) for j in range(n)), "<=",
                            round(act[i] + rng.uniform(0.0, 5.0), 3), name=f"r{i}"))
        else:
            rows.append(Row(tuple((j, A[i][j]) for j in range(n)), ">=",
                            max(0.0, round(act[i] - rng.uniform(0.0, 5.0), 3)), name=f"r{i}"))
    obj = [(j, round(rng.uniform(1.0, 50.0), 2)) for j in range(n)]
    return Model(f"general_ip_{m}x{n}{'' if maximize else '_min'}", vs, rows, obj,
                 "maximize" if maximize else "minimize",
                 meta={"family": "general_ip", "m": m, "n": n, "maximize": maximize,
                       "seed": seed, "planted": [float(v) for v in xh]})


def _gf2_solve(eqs: List[int], rhs: List[int], nvars: int) -> "List[int] | None":
    """Solve a GF(2) system (eqs[i] = bitmask of the variables in equation i).

    Returns one solution (free variables at 0) or None when inconsistent.
    """
    rows = [(e, r) for e, r in zip(eqs, rhs)]
    pivots: List[Tuple[int, int, int]] = []  # (column, mask, rhs)
    for e, r in rows:
        for c, pm, pr in pivots:
            if e >> c & 1:
                e ^= pm
                r ^= pr
        if e == 0:
            if r:
                return None
            continue
        c = e.bit_length() - 1
        # keep the pivot list reduced: eliminate c from earlier pivots
        pivots = [(pc, pm ^ e if pm >> c & 1 else pm, pr ^ r if pm >> c & 1 else pr)
                  for pc, pm, pr in pivots]
        pivots.append((c, e, r))
    x = [0] * nvars
    for c, pm, pr in pivots:
        # reduced form: pm = bit c + free columns (all free columns are 0)
        x[c] = pr
    return x


def lights_out(n: int = 6, infeasible: bool = False, seed: int = 42) -> Model:
    """Lights Out on an n x n board as an integer program.

        min sum x_ij
        s.t. sum_{(u,v) in N[i,j]} x_uv - 2 k_ij = s_ij   (one row per cell)
             x binary (press), k >= 0 integer (no upper bound)
    N[i,j] is the cell and its orthogonal neighbours; s is the initial light
    pattern, and the rows say every light is toggled an odd/even number of
    times as needed to switch it off. Anderson & Feil 1998, "Turning Lights
    Out with Linear Algebra", Math. Mag. 71: solvability is a GF(2) linear
    system. This is the MIPLIB 2017 enlight* class (all rows general linear,
    half binaries, half unbounded integers, objective on the binaries only).
    Feasible mode: s = A p mod 2 for a random press pattern p, so p (and
    k = (A p - s) / 2) is a planted witness; the optimum can be lower.
    infeasible=True: s is drawn outside the GF(2) column space of A, which
    exists only when A is singular (n = 4, 5, 9, 11, ...; Anderson & Feil
    give the nullities). The LP relaxation stays feasible (k continuous
    absorbs any parity), so infeasibility is purely integral. Raises
    ValueError for a nonsingular board.
    """
    rng = random.Random(seed)
    N = n * n
    cell = lambda i, j: i * n + j  # noqa: E731
    nbr = []
    for i in range(n):
        for j in range(n):
            nb = [cell(i, j)]
            for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                if 0 <= i + di < n and 0 <= j + dj < n:
                    nb.append(cell(i + di, j + dj))
            nbr.append(sorted(nb))
    masks = [sum(1 << u for u in nb) for nb in nbr]
    if infeasible:
        for _ in range(1000):
            s = [rng.randint(0, 1) for _ in range(N)]
            if _gf2_solve(masks, s, N) is None:
                break
        else:
            raise ValueError(f"lights_out board {n}x{n} is nonsingular: every pattern solvable")
        press = None
    else:
        press = [rng.randint(0, 1) for _ in range(N)]
        s = [sum(press[u] for u in nb) % 2 for nb in nbr]
    vs = ([Var(0, 1, VTYPE_BIN, name=f"x{c}") for c in range(N)]
          + [Var(0.0, INF, VTYPE_INT, name=f"k{c}") for c in range(N)])
    rows = [Row(tuple([(u, 1.0) for u in nbr[c]] + [(N + c, -2.0)]), "==", float(s[c]),
                name=f"cell{c}") for c in range(N)]
    obj = [(c, 1.0) for c in range(N)]
    meta = {"family": "lights_out", "n": n, "infeasible": infeasible, "seed": seed}
    if press is not None:
        k = [(sum(press[u] for u in nbr[c]) - s[c]) // 2 for c in range(N)]
        meta["planted"] = [float(v) for v in press] + [float(v) for v in k]
    return Model(f"lights_out_{n}{'_infeas' if infeasible else ''}", vs, rows, obj, "minimize",
                 meta=meta)


def mixed_int_knapsack(m: int = 12, n_int: int = 24, n_bin: int = 10, seed: int = 42) -> Model:
    """Multi-row mixed-integer knapsack with bounded integers, binaries and slacks.

        max  sum p_j x_j + sum q_k z_k - sum h_i s_i
        s.t. sum_j a_ij x_j + sum_k g_ik z_k - s_i <= b_i   (general linear)
             x_j <= u'_j                                     (singleton rows)
             0 <= x_j <= u_j integer, z binary, s >= 0
    Atamtürk 2003, "On the facets of the mixed-integer knapsack polyhedron",
    Math. Prog. 98: the mixed-integer knapsack set (bounded general integers
    plus one continuous variable per row) whose test sets are the MIPLIB
    2017 mik.* instances. The singleton rows restate tighter integer upper
    bounds as rows, as those instances do (mik-250-20-75-4: 61% singleton
    rows); presolve should turn them into bounds. Coefficients in [1, 2000],
    rhs at 40% of the row's bounded maximum, so slacks and branching matter.
    Witness: everything at 0 (b >= 0).
    """
    rng = random.Random(seed)
    nx, nz = n_int, n_bin
    u = [rng.randint(5, 20) for _ in range(nx)]
    ut = [rng.randint(2, uj - 1) for uj in u]
    xi = lambda j: j  # noqa: E731
    zi = lambda k: nx + k  # noqa: E731
    si = lambda i: nx + nz + i  # noqa: E731
    vs = ([Var(0.0, float(u[j]), VTYPE_INT, name=f"x{j}") for j in range(nx)]
          + [Var(0, 1, VTYPE_BIN, name=f"z{k}") for k in range(nz)]
          + [Var(0.0, INF, name=f"s{i}") for i in range(m)])
    rows = []
    for i in range(m):
        cols_x = rng.sample(range(nx), max(2, nx // 2))
        cols_z = rng.sample(range(nz), max(1, nz // 2)) if nz else []
        coeffs = [(xi(j), float(rng.randint(1, 2000))) for j in sorted(cols_x)]
        coeffs += [(zi(k), float(rng.randint(1, 2000))) for k in sorted(cols_z)]
        top = sum(c * (ut[j] if j < nx else 1) for j, c in coeffs)
        coeffs.append((si(i), -1.0))
        rows.append(Row(tuple(coeffs), "<=", float(int(0.4 * top)), name=f"knap{i}"))
    rows += [Row(((xi(j), 1.0),), "<=", float(ut[j]), name=f"ub{j}") for j in range(nx)]
    obj = ([(xi(j), float(rng.randint(100, 2000))) for j in range(nx)]
           + [(zi(k), float(rng.randint(100, 4000))) for k in range(nz)]
           + [(si(i), float(rng.randint(1, 5))) for i in range(m)])
    obj = [(j, c if j < nx + nz else -c) for j, c in obj]
    return Model(f"mixed_int_knapsack_{m}x{nx}i{nz}b", vs, rows, obj, "maximize",
                 meta={"family": "mixed_int_knapsack", "m": m, "n_int": nx, "n_bin": nz,
                       "seed": seed, "planted": [0.0] * len(vs)})


def p_median(n: int = 30, p: int = 5, bounded: bool = False, seed: int = 42) -> Model:
    """p-median, Elloumi's radius formulation (mostly continuous, few binaries).

        min  sum_i sum_{k>=2} (D_i^k - D_i^{k-1}) z_ik
        s.t. z_ik + sum_{j : d_ij < D_i^k} y_j >= 1   (i, k >= 2)
             sum_j y_j = p                          (cardinality)
             z >= 0, y binary
    D_i^1 = 0 < D_i^2 < ... are the distinct distances from point i; z_ik = 1
    when no open median lies closer than D_i^k. Elloumi 2010, "A tighter
    formulation of the p-median problem", J. Comb. Optim. 19 (after
    Cornuéjols, Nemhauser & Wolsey 1980, "A canonical representation of
    simple plant location problems and its applications", SIAM J. Algebraic
    Discrete Methods 1). This is the MIPLIB 2017 rmatr* shape: 99% continuous
    unbounded columns, a few binaries, rows mixing one continuous with many
    binaries. Points uniform in [0, 100]^2, integral Euclidean distances.
    bounded=True states the implied z_ik <= 1 as a column bound (boxed
    continuous columns, the istanbul/pg shape; same optimum).
    Witness: a random set of p medians with its exact cost (planted_obj).
    """
    rng = random.Random(seed)
    pts = [(rng.uniform(0, 100), rng.uniform(0, 100)) for _ in range(n)]
    d = [[int(round(((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5)) for b in pts] for a in pts]
    vs = [Var(0, 1, VTYPE_BIN, name=f"y{j}") for j in range(n)]
    rows, obj = [], []
    zcols = []  # (i, D_prev, D_k, col)
    for i in range(n):
        D = sorted(set(d[i]))
        for k in range(1, len(D)):
            c = len(vs)
            vs.append(Var(0.0, 1.0 if bounded else INF, name=f"z{i}_{k}"))
            zcols.append((i, D[k], c))
            closer = [j for j in range(n) if d[i][j] < D[k]]
            rows.append(Row(tuple([(c, 1.0)] + [(j, 1.0) for j in closer]), ">=", 1.0,
                            name=f"rad{i}_{k}"))
            obj.append((c, float(D[k] - D[k - 1])))
    rows.append(Row(tuple((j, 1.0) for j in range(n)), "==", float(p), name="card"))
    med = set(rng.sample(range(n), p))
    planted = [1.0 if j in med else 0.0 for j in range(n)]
    planted += [0.0] * (len(vs) - n)
    for i, Dk, c in zcols:
        if not any(d[i][j] < Dk for j in med):
            planted[c] = 1.0
    pobj = float(sum(min(d[i][j] for j in med) for i in range(n)))
    return Model(f"p_median_{n}p{p}{'_ub' if bounded else ''}", vs, rows, obj, "minimize",
                 meta={"family": "p_median", "n": n, "p": p, "bounded": bounded, "seed": seed,
                       "planted": planted, "planted_obj": pobj})


def elastic_cover(n: int = 40, m: int = 6, card_frac: float = 0.3, seed: int = 42) -> Model:
    """Multi-row binary covering knapsack with one elastic (penalised) shortfall.

        min  sum c_j x_j + P s
        s.t. sum_j a_ij x_j + s >= b_i     (mixed binary, one shared continuous)
             sum_j x_j <= K                (invariant knapsack)
             x binary, s >= 0
    The elastic row form of goal programming (Charnes & Cooper 1961,
    "Management Models and Industrial Applications of Linear Programming",
    Wiley): one continuous variable absorbs the shortfall of every row at a
    penalty P. This is the MIPLIB 2017 mas74/mas76 shape (12-13 rows x 150
    binaries + 1 continuous, dense rows with coefficients up to 1e4, every
    column in the objective, objective dynamism ~1e5). a_ij ~ U{1..9999}, c_j ~ U{1..2000},
    b_i = 0.5 sum_j a_ij, K = card_frac n, P = 1e5.
    Witness: x = 0, s = max_i b_i.
    """
    rng = random.Random(seed)
    a = [[rng.randint(1, 9999) for _ in range(n)] for _ in range(m)]
    b = [int(0.5 * sum(r)) for r in a]
    K = max(1, int(card_frac * n))
    vs = _bins(n) + [Var(0.0, INF, name="s")]
    rows = [Row(tuple([(j, float(a[i][j])) for j in range(n)] + [(n, 1.0)]), ">=", float(b[i]),
                name=f"cov{i}") for i in range(m)]
    rows.append(Row(tuple((j, 1.0) for j in range(n)), "<=", float(K), name="card"))
    obj = [(j, float(rng.randint(1, 2000))) for j in range(n)] + [(n, 1e5)]
    planted = [0.0] * n + [float(max(b))]
    return Model(f"elastic_cover_{n}x{m}", vs, rows, obj, "minimize",
                 meta={"family": "elastic_cover", "n": n, "m": m, "K": K, "seed": seed,
                       "planted": planted})


# Per-block sizes for block_angular: "tiny" keeps a 2-block instance inside the
# mip_diff small-case budget, "small" is about MIP_SMALL, "block" is the
# default/medium block.
BLOCK_SIZES = {
    "tiny": {"gap": dict(n_agents=2, n_jobs=4), "lot_sizing": dict(T=4),
             "multi_knapsack": dict(n=8, m=2), "fixed_charge_network": dict(nodes=4, arcs=6),
             "facility_mixed": dict(n_fac=3, n_cust=3, hubs=1)},
    "small": {"gap": dict(n_agents=3, n_jobs=8), "lot_sizing": dict(T=8),
              "multi_knapsack": dict(n=20, m=5), "fixed_charge_network": dict(nodes=6, arcs=12),
              "facility_mixed": dict(n_fac=4, n_cust=5, hubs=1)},
    "block": {"gap": dict(n_agents=3, n_jobs=12), "lot_sizing": dict(T=15),
              "multi_knapsack": dict(n=25, m=3), "fixed_charge_network": dict(nodes=8, arcs=20),
              "facility_mixed": dict(n_fac=5, n_cust=8, hubs=1)},
}
BLOCK_MIX = ("gap", "lot_sizing", "multi_knapsack", "fixed_charge_network", "facility_mixed")


def block_angular(blocks: int = 4, base: str = "mixed", block_size: str = "block",
                  link_rows: int = 3, link_vars: int = 4, infeasible: bool = False,
                  seed: int = 42) -> Model:
    """Block-angular MIP: independent sub-models coupled by a few linking rows.

        min  sum_b c_b' x_b
        s.t. A_b x_b (sense) b_b          for each block b  (the base family's rows)
             sum_b D_b x_b (<= or >=) e   (link_rows coupling rows)
    The block-angular shape of Dantzig & Wolfe 1960, "Decomposition principle
    for linear programs", Oper. Res. 8; detecting it in MIPLIB instances is
    the subject of Bergner et al. 2015, "Automatic Dantzig-Wolfe reformulation
    of mixed integer programs", Math. Prog. 149. MIPLIB 2017 has instances
    with several components (neos-3004026-krka: 2146) that no other family
    produced.
    base="mixed" cycles the blocks through BLOCK_MIX, so one model holds
    assignment, lot-sizing, knapsack, fixed-charge-flow and facility rows;
    otherwise every block is `base`. Maximize blocks enter with a negated
    objective. Linking rows pick link_vars columns per block, with coefficients
    in [1, 20]. They alternate <= (shared resource: rhs = planted activity +
    10% + 1) and >= (minimum joint activity: rhs = 0.9 x planted activity), so
    the planted witness (the blocks' witnesses side by side) stays feasible.
    link_rows=0 gives `blocks` independent components.
    infeasible=True replaces the last block with gap(infeasible=True).
    That block has a feasible LP but no integer point, so the model is
    infeasible while the other blocks and the links stay LP-feasible.
    """
    rng = random.Random(seed)
    sizes = BLOCK_SIZES[block_size]
    names = [BLOCK_MIX[b % len(BLOCK_MIX)] if base == "mixed" else base for b in range(blocks)]
    vs: List[Var] = []
    rows: List[Row] = []
    obj: List[Tuple[int, float]] = []
    planted: List[float] = []
    feasible_cols: List[List[int]] = []
    for b, fam in enumerate(names):
        kw = dict(sizes.get(fam, MIP_SMALL.get(fam, {})))
        kw["seed"] = rng.randrange(1 << 30)
        bad = infeasible and b == blocks - 1
        if bad:
            sub = gap(**{**sizes["gap"], "seed": kw["seed"], "infeasible": True})
        else:
            sub = REGISTRY[fam](**kw)
        off = len(vs)
        sign = -1.0 if sub.obj_sense == "maximize" else 1.0
        vs += [Var(v.lb, v.ub, v.vtype, name=f"b{b}_{v.name}") for v in sub.vars]
        rows += [Row(tuple((off + j, c) for j, c in r.coeffs), r.sense, r.rhs,
                     name=f"b{b}_{r.name}") for r in sub.rows]
        obj += [(off + j, sign * c) for j, c in sub.obj]
        p = sub.meta.get("planted")
        x = [0.0] * sub.ncols
        if p is not None:
            for j, v in (p.items() if isinstance(p, dict) else enumerate(p)):
                x[int(j)] = float(v)
        planted += x
        if not bad:
            feasible_cols.append(list(range(off, off + sub.ncols)))
    for r in range(link_rows):
        coeffs = []
        for cols in feasible_cols:
            for j in rng.sample(cols, min(link_vars, len(cols))):
                coeffs.append((j, float(rng.randint(1, 20))))
        act = sum(c * planted[j] for j, c in coeffs)
        if r % 2 == 0:
            rows.append(Row(tuple(coeffs), "<=", float(int(1.1 * act) + 1), name=f"link{r}"))
        else:
            rows.append(Row(tuple(coeffs), ">=", float(int(0.9 * act)), name=f"link{r}"))
    meta = {"family": "block_angular", "blocks": blocks, "base": base, "block_size": block_size,
            "link_rows": link_rows, "seed": seed, "infeasible": infeasible,
            "block_families": names}
    if not infeasible:
        meta["planted"] = planted
    tag = "_infeas" if infeasible else ""
    return Model(f"block_angular_{blocks}x{base}_{block_size}_l{link_rows}{tag}", vs, rows, obj,
                 "minimize", meta=meta)


REGISTRY = {
    "degenerate_ladder": degenerate_ladder,
    "set_cover": set_cover,
    "network_flow": network_flow,
    "knapsack_cover": knapsack_cover,
    "corner_lp": corner_lp,
    "transport_degen": transport_degen,
    "set_packing": set_packing,
    "set_partition": set_partition,
    "multi_knapsack": multi_knapsack,
    "gen_int_knapsack": gen_int_knapsack,
    "cap_facility_location": cap_facility_location,
    "uncap_facility_location": uncap_facility_location,
    "lot_sizing": lot_sizing,
    "gap": gap,
    "fixed_charge_network": fixed_charge_network,
    "bin_packing": bin_packing,
    "graph_coloring": graph_coloring,
    "max_indep_set": max_indep_set,
    "market_split": market_split,
    "jobshop_bigM": jobshop_bigM,
    "network_design": network_design,
    "facility_mixed": facility_mixed,
    "budget_max_cover": budget_max_cover,
    "int_production": int_production,
    "block_angular": block_angular,
    "general_ip": general_ip,
    "lights_out": lights_out,
    "mixed_int_knapsack": mixed_int_knapsack,
    "p_median": p_median,
    "elastic_cover": elastic_cover,
}

# First positional arg that --size maps to, per family.
SIZE_ARG = {
    "degenerate_ladder": "n",
    "set_cover": "m",
    "network_flow": "nodes",
    "knapsack_cover": "n",
    "corner_lp": "n",
    "transport_degen": "m",
    "set_packing": "m",
    "set_partition": "m",
    "multi_knapsack": "n",
    "gen_int_knapsack": "n",
    "cap_facility_location": "n_cust",
    "uncap_facility_location": "n_cust",
    "lot_sizing": "T",
    "gap": "n_jobs",
    "fixed_charge_network": "nodes",
    "bin_packing": "n_items",
    "graph_coloring": "n_nodes",
    "max_indep_set": "n_nodes",
    "market_split": "m",
    "jobshop_bigM": "n_jobs",
    "network_design": "nodes",
    "facility_mixed": "n_cust",
    "budget_max_cover": "n_sets",
    "int_production": "periods",
    "block_angular": "blocks",
    "general_ip": "n",
    "lights_out": "n",
    "mixed_int_knapsack": "n_int",
    "p_median": "n",
    "elastic_cover": "n",
}

# MIP families (mip-plan §2a) in registry order.
MIP_FAMILIES = (
    "set_cover", "knapsack_cover", "set_packing", "set_partition", "multi_knapsack",
    "gen_int_knapsack", "cap_facility_location", "uncap_facility_location", "lot_sizing",
    "gap", "fixed_charge_network", "bin_packing", "graph_coloring", "max_indep_set",
    "market_split", "jobshop_bigM",
    # mixed-structure families
    "network_design", "facility_mixed", "budget_max_cover", "int_production", "block_angular",
    "general_ip", "lights_out", "mixed_int_knapsack", "p_median", "elastic_cover",
)

# `small` size per MIP family: <= ~30 vars, <= ~25 rows, so a slow exact B&B
# (and the mip_diff gate) can solve every instance. Variant kwargs (form,
# equality, infeasible) are layered on top by the caller.
MIP_SMALL = {
    "set_cover": dict(m=12, n=24, density=0.2),
    "knapsack_cover": dict(n=20),
    "set_packing": dict(m=12, n=25, density=0.2),
    "set_partition": dict(m=10, n=28),
    "multi_knapsack": dict(n=20, m=5),
    "gen_int_knapsack": dict(n=12, m=2),
    "cap_facility_location": dict(n_fac=4, n_cust=6, strong=False),
    "uncap_facility_location": dict(n_fac=3, n_cust=6),
    "lot_sizing": dict(T=8),
    "gap": dict(n_agents=3, n_jobs=8),
    "fixed_charge_network": dict(nodes=6, arcs=12),
    "bin_packing": dict(n_items=5, extra_bins=0),
    "graph_coloring": dict(n_nodes=6, density=0.4, max_edges=5),
    "max_indep_set": dict(n_nodes=15, density=0.2),
    "market_split": dict(m=2, n=10),
    "jobshop_bigM": dict(n_jobs=3, n_machines=3),
    "network_design": dict(nodes=4, arcs=7, commodities=2, strong=False),
    "facility_mixed": dict(n_fac=4, n_cust=4, hubs=1),
    "budget_max_cover": dict(n_sets=8, n_elems=12, density=0.25, k=3, conflicts=2),
    "int_production": dict(products=2, periods=3),
    "block_angular": dict(blocks=2, block_size="tiny", link_rows=1, link_vars=3),
    "general_ip": dict(m=4, n=6),
    "lights_out": dict(n=4),
    "mixed_int_knapsack": dict(m=4, n_int=8, n_bin=4),
    "p_median": dict(n=5, p=2),
    "elastic_cover": dict(n=16, m=3),
}
