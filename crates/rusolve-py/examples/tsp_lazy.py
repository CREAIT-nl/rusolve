"""Lazy-constraint TSP in Python: `python crates/rusolve-py/examples/tsp_lazy.py 9 4`.

The model has only the degree-2 rows; the callback adds a subtour-elimination row for every
subtour of each candidate (Dantzig, Fulkerson, Johnson, Operations Research 2 (1954) 393-410).
Mirrors `crates/rusolve-ekk/examples/tsp_lazy.rs`.
"""
import sys

import hooks_common as hc

n = int(sys.argv[1]) if len(sys.argv) > 1 else 9
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 4
pts = hc.cities(n, seed)
m = hc.tsp_model(pts)
sol = m.solve(callback=hc.SubtourElimination(n), threads=1)
stats = m.callback_stats
assert sol.x is not None and stats is not None
print(f"n={n} seed={seed} status={sol.status} tour={sol.objective:.6f} brute={hc.brute_force_tour(pts):.6f}")
print(f"rounds={stats['lazy']['rounds']} rows_added={stats['lazy']['rows_added']} model rows={m.num_constrs}")
