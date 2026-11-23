"""Prototype a branching rule and a node selector in Python:
`python crates/rusolve-py/examples/prototype_branching.py`.

Both hooks can only change the order of the search, never the optimum; the tree size changes.
Run per-class components here first, then promote the winner to Rust or to a Wasm component.
"""
import hooks_common as hc
import rusolve as rs


class HighIndexFirst:
    def branching_score(self, ctx: rs.CallbackContext) -> list[float]:
        assert ctx.candidates is not None
        return [float(col) for col, _value in ctx.candidates]

    def select_node(self, ctx: rs.CallbackContext) -> int | None:
        nodes = ctx.open_nodes
        assert nodes is not None
        return max(range(len(nodes)), key=lambda k: nodes[k][2])


base = hc.multi_knapsack(34, 5)
b = base.solve(presolve=False, threads=1, options=hc.TREE_OPTIONS)
m = hc.multi_knapsack(34, 5)
h = m.solve(callback=HighIndexFirst(), threads=1, options=hc.TREE_OPTIONS)
print(f"default  objective {b.objective} nodes {b.nodes}")
print(f"python   objective {h.objective} nodes {h.nodes}  stats {m.callback_stats}")
