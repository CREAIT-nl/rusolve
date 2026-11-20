"""Instances shared by the callback examples and `tests/test_callbacks.py`.

The generators reproduce the Rust ones in `crates/rusolve-ekk/examples/common/mod.rs`
(same LCG, same seeds), so the Python callbacks solve exactly the instances that
`crates/rusolve-ekk/tests/lazy_cuts.rs` checks.
"""
from __future__ import annotations

import itertools
import math
from typing import Any

import rusolve as rs

MASK = (1 << 64) - 1


class Lcg:
    def __init__(self, seed: int) -> None:
        self.s = seed & MASK

    def next(self) -> int:
        self.s = (self.s * 6364136223846793005 + 1442695040888963407) & MASK
        return self.s >> 33


# ---------------------------------------------------------------- TSP (lazy constraints)


def cities(n: int, seed: int) -> list[tuple[float, float]]:
    g = Lcg(seed * 7919 + 17)
    pts = []
    for _ in range(n):
        x = float(g.next() % 100)
        y = float(g.next() % 100)
        pts.append((x, y))
    return pts


def dist(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def edge_col(n: int, i: int, j: int) -> int:
    """Column of the undirected edge {i, j} (edges ordered (0,1) (0,2) ... (1,2) ...)."""
    if i > j:
        i, j = j, i
    return i * n - i * (i + 1) // 2 + (j - i - 1)


def tsp_model(pts: list[tuple[float, float]]) -> rs.Model:
    """Symmetric TSP with only the degree-2 rows: subtours are the callback's job."""
    n = len(pts)
    m = rs.Model(name=f"tsp{n}")
    x = m.add_vars(n * (n - 1) // 2, binary=True)
    m.minimize(
        rs.quicksum(dist(pts[i], pts[j]) * x[edge_col(n, i, j)] for i in range(n) for j in range(i + 1, n))
    )
    for i in range(n):
        m.add_constr(rs.quicksum(x[edge_col(n, i, j)] for j in range(n) if j != i) == 2)
    return m


def components(sol: list[float], n: int) -> list[list[int]]:
    """Connected components of the chosen edges (value > 0.5)."""
    seen = [False] * n
    out: list[list[int]] = []
    for s in range(n):
        if seen[s]:
            continue
        seen[s] = True
        comp, stack = [], [s]
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in range(n):
                if v != u and not seen[v] and sol[edge_col(n, u, v)] > 0.5:
                    seen[v] = True
                    stack.append(v)
        out.append(sorted(comp))
    return out


class SubtourElimination:
    """Lazy constraints: one subtour-elimination row per subtour of the candidate."""

    def __init__(self, n: int) -> None:
        self.n = n

    def lazy_constraints(self, ctx: rs.CallbackContext) -> list[tuple[list[int], list[float], float]]:
        assert ctx.solution is not None
        comps = components(ctx.solution, self.n)
        if len(comps) == 1:
            return []
        cuts = []
        for s in comps:
            if len(s) >= 3:
                cols = [edge_col(self.n, i, j) for a, i in enumerate(s) for j in s[a + 1 :]]
                cuts.append((cols, [1.0] * len(cols), len(s) - 1.0))
        return cuts


def brute_force_tour(pts: list[tuple[float, float]]) -> float:
    """Optimal tour length by Held-Karp dynamic programming (exact, any n up to about 15)."""
    n = len(pts)
    d = [[dist(pts[i], pts[j]) for j in range(n)] for i in range(n)]
    full = 1 << (n - 1)
    inf = float("inf")
    best = [[inf] * (n - 1) for _ in range(full)]
    for j in range(n - 1):
        best[1 << j][j] = d[0][j + 1]
    for mask in range(full):
        for j in range(n - 1):
            c = best[mask][j]
            if c == inf or not mask >> j & 1:
                continue
            for k in range(n - 1):
                if not mask >> k & 1:
                    nm = mask | 1 << k
                    v = c + d[j + 1][k + 1]
                    if v < best[nm][k]:
                        best[nm][k] = v
    return min(best[full - 1][j] + d[j + 1][0] for j in range(n - 1))


def brute_force_permutations(pts: list[tuple[float, float]]) -> float:
    """Same optimum by trying every permutation (city 0 fixed); for n <= 8 in tests."""
    n = len(pts)
    best = float("inf")
    for p in itertools.permutations(range(1, n)):
        tour = (0, *p, 0)
        best = min(best, sum(dist(pts[a], pts[b]) for a, b in zip(tour, tour[1:])))
    return best


# ---------------------------------------------------------------- knapsack covers (user cuts)


def knapsack_data() -> tuple[list[float], list[float], float]:
    g = Lcg(4242)
    n = 14
    weights = [float(10 + g.next() % 30) for _ in range(n)]
    values = [w + float(g.next() % 7) for w in weights]
    return values, weights, math.floor(sum(weights) * 0.4)


def minimal_covers(weights: list[float], cap: float) -> list[list[int]]:
    """Every minimal cover C (sum_C w > cap, and dropping any item fits)."""
    n = len(weights)
    out = []
    for mask in range(1, 1 << n):
        items = [i for i in range(n) if mask >> i & 1]
        total = sum(weights[i] for i in items)
        if total > cap and all(total - weights[i] <= cap for i in items):
            out.append(items)
    return out


def knapsack_model(covers_in_model: list[list[int]] | None = None) -> rs.Model:
    values, weights, cap = knapsack_data()
    m = rs.Model(sense="max", name="knapsack")
    x = m.add_vars(len(values), binary=True)
    m.add_constr(rs.quicksum(w * x[i] for i, w in enumerate(weights)) <= cap)
    for cv in covers_in_model or []:
        m.add_constr(rs.quicksum(x[i] for i in cv) <= len(cv) - 1)
    m.maximize(rs.quicksum(v * x[i] for i, v in enumerate(values)))
    return m


class CoverSeparator:
    """User cuts: the minimal cover inequalities `sum_C x <= |C| - 1`."""

    def __init__(self, covers: list[list[int]]) -> None:
        self.covers = covers

    def separate(self, ctx: rs.CallbackContext) -> list[tuple[list[int], list[float], float]]:
        return [(cv, [1.0] * len(cv), len(cv) - 1.0) for cv in self.covers]


# ---------------------------------------------------------------- a tree-forcing knapsack


def multi_knapsack(n: int = 28, rows: int = 4) -> rs.Model:
    """The 28-item 4-row 0/1 knapsack of `crates/rusolve-ekk/tests/hooks_mip.rs` (needs a tree)."""
    g = Lcg(12345)
    values = [float(20 + g.next() % 60) for _ in range(n)]
    weights = [[float(10 + g.next() % 50) for _ in range(n)] for _ in range(rows)]
    m = rs.Model(sense="max", name="multiknap")
    x = m.add_vars(n, binary=True)
    for w in weights:
        m.add_constr(rs.quicksum(w[j] * x[j] for j in range(n)) <= math.floor(sum(w) * 0.35))
    m.maximize(rs.quicksum(values[j] * x[j] for j in range(n)))
    return m


TREE_OPTIONS: dict[str, Any] = {"mip_heuristic_effort": 0.0, "mip_heuristic_run_feasibility_jump": False}

# the cover-cut instance solves at the root unless the heuristics are off and the gap is closed
CUT_OPTIONS: dict[str, Any] = {**TREE_OPTIONS, "mip_rel_gap": 0.0, "mip_abs_gap": 0.0}
