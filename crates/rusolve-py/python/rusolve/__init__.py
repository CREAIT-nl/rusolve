"""rusolve: a fast LP/MIP solver with a simple Python model API.

    import rusolve as rs
    m = rs.Model()
    x = m.add_vars(3, ub=1, binary=True)
    m.add_constr(x[0] + x[1] + x[2] <= 2)
    m.maximize(5 * x[0] + 4 * x[1] + 3 * x[2])
    sol = m.solve()
    sol.status, sol.objective, sol.x

See docs/design/python-api.md. `rusolve.Highs` is the highspy-shaped compatibility class;
`rusolve.scipy` provides `milp` and `linprog` with SciPy's signatures.
"""
from . import _rusolve
from ._rusolve import *  # noqa: F401,F403
from .errors import *  # noqa: F401,F403
from . import errors

inf = float("inf")

__all__ = [n for n in dir(_rusolve) if not n.startswith("_")] + list(errors.__all__) + ["errors", "inf"]
