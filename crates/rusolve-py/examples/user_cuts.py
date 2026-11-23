"""User cuts in Python: `python crates/rusolve-py/examples/user_cuts.py`.

The minimal cover inequalities of a 14-item knapsack, as model rows and from a Python
`separate` callback, reach the plain optimum (180). Mirrors `crates/rusolve-ekk/examples/user_cuts.rs`.
"""
import hooks_common as hc

values, weights, cap = hc.knapsack_data()
covers = hc.minimal_covers(weights, cap)
print(len(covers), "minimal covers")
print("plain          ", hc.knapsack_model().solve(presolve=False, threads=1, options=hc.CUT_OPTIONS).objective)
print("cuts in model  ", hc.knapsack_model(covers).solve(presolve=False, threads=1, options=hc.CUT_OPTIONS).objective)
m = hc.knapsack_model()
s = m.solve(callback=hc.CoverSeparator(covers), threads=1, options=hc.CUT_OPTIONS)
stats = m.callback_stats
assert stats is not None
print("python separator", s.objective, "accepted", stats["cuts_accepted"], "cuts")
