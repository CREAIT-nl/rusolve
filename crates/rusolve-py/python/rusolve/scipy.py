"""SciPy-compatible ``milp`` and ``linprog`` on top of ``rusolve.Model``.

Switching from SciPy is a change of import:

    from rusolve.scipy import milp, linprog

Signatures, defaults, result fields and status codes follow ``scipy.optimize`` (0 optimal,
1 iteration/time/node limit, 2 infeasible, 3 unbounded, 4 other). ``Bounds`` and
``LinearConstraint`` are duck-typed (``.lb``, ``.ub``, ``.A``), so scipy itself is never
imported here. Matrices are handed to ``Model.add_matrix`` as given (dense or scipy.sparse):
there is no densify step. See docs/design/python-api.md, "Compatibility layers".
"""
import warnings

import numpy as np

from . import Model

_INF = float("inf")

__all__ = ["OptimizeResult", "linprog", "milp"]


class OptimizeResult(dict):
    """A dict with attribute access, like ``scipy.optimize.OptimizeResult``."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if not self:
            return self.__class__.__name__ + "()"
        w = max(map(len, self)) + 1
        return "\n".join(k.rjust(w) + ": " + repr(v) for k, v in sorted(self.items()))

    def __dir__(self):
        return list(self.keys())


# status string (closed result-status enum) -> (scipy status, scipy message)
_LIMITS = {"TimeLimit": "Time limit reached. ", "IterationLimit": "Iteration limit reached. ",
           "NodeLimit": "Node limit reached. ", "SolutionLimit": "Solution limit reached. "}


def _status(s):
    if s == "Optimal":
        return 0, "Optimization terminated successfully. "
    if s in _LIMITS:
        return 1, _LIMITS[s]
    if s == "Infeasible":
        return 2, "The problem is infeasible. "
    if s == "Unbounded":
        return 3, "The problem is unbounded. "
    if s == "InfeasibleOrUnbounded":
        return 4, "The problem is unbounded or infeasible. "
    return 4, "The solver did not provide a recognised status. "


# ----------------------------------------------------------------------------- shared helpers
_MILP_OPTIONS = {"disp", "presolve", "time_limit", "node_limit", "mip_rel_gap"}


def _solve_options(options, supported, ipm=False):
    """SciPy option dict -> (kwargs for Model.solve, raw HiGHS options)."""
    kw, raw = {}, {}
    if ipm:
        raw["solver"] = "ipm"
    if not options:
        return kw, raw
    options = dict(options)
    unknown = set(options) - supported
    if unknown:
        warnings.warn(f"Unrecognized options detected: {unknown}. These will be passed to "
                      "the solver verbatim.", RuntimeWarning, stacklevel=3)
    if "disp" in options:
        kw["verbose"] = bool(options.pop("disp"))
    if "time_limit" in options:
        v = options.pop("time_limit")
        if v is not None:
            kw["time_limit"] = float(v)
    if "presolve" in options:
        v = options.pop("presolve")
        if v is not None:
            kw["presolve"] = bool(v)
    if "mip_rel_gap" in options:
        v = options.pop("mip_rel_gap")
        if v is not None:
            kw["mip_gap"] = float(v)
    if "node_limit" in options:
        v = options.pop("node_limit")
        if v is not None:
            raw["mip_max_nodes"] = int(v)
    if "maxiter" in options:
        v = options.pop("maxiter")
        if v is not None:
            raw["simplex_iteration_limit"] = int(v)
            raw["ipm_iteration_limit"] = int(v)
    raw.update(options)
    return kw, raw


def _is_sparse(a):
    return hasattr(a, "indptr") or hasattr(a, "tocsr")


def _as_matrix(A):
    """Sparse or dense array as given; lists become float64 arrays."""
    if _is_sparse(A):
        return A
    if isinstance(A, np.ndarray) and A.dtype == np.float64:
        return A if A.ndim == 2 else np.atleast_2d(A)
    return np.atleast_2d(np.asarray(A, dtype=np.float64))


def _vec(v, n):
    return np.broadcast_to(np.asarray(v, dtype=np.float64), (n,))


def _add_rows(m, A, lo, hi):
    """Add ``lo <= A x <= hi`` rows (one ranged add_matrix call; infinite sides are open)."""
    if A.shape[0]:
        m.add_matrix(A, lo=lo, hi=hi)


# ----------------------------------------------------------------------------------- milp
def _constraint_parts(constraints):
    """-> list of (A, lb, ub) from a LinearConstraint, a sequence of them, or (A, lb, ub) tuples."""
    if not constraints and not (hasattr(constraints, "shape") and constraints.size):
        return []
    if hasattr(constraints, "A"):
        constraints = [constraints]
    out = []
    for con in constraints:
        if hasattr(con, "A"):
            A, lb, ub = con.A, con.lb, con.ub
        else:
            A, lb, ub = con
        A = _as_matrix(A)
        k = A.shape[0]
        out.append((A, _vec(lb, k), _vec(ub, k)))
    return out


def _bnd(v, n):
    """Scalar bound as is; array bound broadcast to length n (Bounds(0, 1).lb has shape (1,))."""
    if np.ndim(v) == 0:
        return v
    a = np.asarray(v, dtype=np.float64)
    return a if a.shape == (n,) else np.broadcast_to(a, (n,))


def _bounds_milp(bounds, n):
    if bounds is None:
        return 0.0, _INF
    if hasattr(bounds, "lb"):
        return _bnd(bounds.lb, n), _bnd(bounds.ub, n)
    lb, ub = bounds  # SciPy does Bounds(*bounds)
    return _bnd(lb, n), _bnd(ub, n)


def _integrality(integrality, n):
    if integrality is None:
        return False
    a = np.asarray(integrality)
    if a.size and (a.min() < 0 or a.max() > 3):
        raise ValueError("`integrality` must contain integers 0-3 and be broadcastable "
                         "to `c.shape`.")
    if a.ndim == 0:
        return bool(a == 1 or a == 3)
    a = np.broadcast_to(a, (n,))
    integ = (a == 1) | (a == 3)
    if not integ.any():
        return False
    if integ.all():
        return True
    return integ.tolist()


def _semi(integrality, n, m):
    """Mark the semi-continuous (2) and semi-integer (3) columns of ``integrality`` on ``m``."""
    if integrality is None:
        return
    a = np.broadcast_to(np.asarray(integrality), (n,))
    for code, integer in ((2, False), (3, True)):
        idx = np.flatnonzero(a == code)
        if idx.size:
            m.set_semi(idx.tolist(), integer=integer)


def _check_vec(c):
    try:
        n = len(c)
    except TypeError:
        c = [c]
        n = 1
    if n == 0:
        raise ValueError("`c` must be a one-dimensional array of finite numbers "
                         "with at least one element.")
    return c, n


def milp(c, *, integrality=None, bounds=None, constraints=None, options=None):
    """Mixed-integer linear programming with the signature of ``scipy.optimize.milp``.

    minimize ``c @ x`` subject to ``lb <= A @ x <= ub`` constraints, variable bounds
    (default ``0 <= x < inf``) and ``integrality`` (0 continuous, 1 integer, 2 semi-continuous, 3 semi-integer).
    Returns an ``OptimizeResult`` with ``x``, ``fun``, ``status``, ``success``, ``message``,
    ``mip_node_count``, ``mip_dual_bound``, ``mip_gap``.
    """
    c, n = _check_vec(c)
    kw, raw = _solve_options(options, _MILP_OPTIONS)
    m = Model()
    lb, ub = _bounds_milp(bounds, n)
    m.add_vars(n, lb=lb, ub=ub, cost=c, integer=_integrality(integrality, n))
    _semi(integrality, n, m)
    for A, lo, hi in _constraint_parts(constraints):
        if A.shape[1] != n:
            raise ValueError("The shape of `A` must be (len(b_l), len(c)).")
        _add_rows(m, A, lo, hi)
    sol = m.solve(options=raw or None, **kw)
    status, msg = _status(sol.status)
    x = sol.x if sol.ok else None
    return OptimizeResult(
        x=None if x is None else np.asarray(x, dtype=np.float64),
        fun=sol.objective if x is not None else None,
        status=status, message=f"{msg}(rusolve status: {sol.status})", success=status == 0,
        mip_node_count=sol.nodes, mip_dual_bound=sol.bound, mip_gap=sol.gap)


# -------------------------------------------------------------------------------- linprog
_LINPROG_OPTIONS = {"disp", "presolve", "time_limit", "maxiter", "dual_feasibility_tolerance",
                    "primal_feasibility_tolerance", "ipm_optimality_tolerance", "mip_rel_gap",
                    "node_limit", "simplex_dual_edge_weight_strategy", "unknown_options"}
_METHODS = {"highs": False, "highs-ds": False, "highs-ipm": True}


def _linprog_bounds(bounds, n):
    if bounds is None:
        return 0.0, _INF
    if hasattr(bounds, "lb"):
        return _bnd(bounds.lb, n), _bnd(bounds.ub, n)
    b = list(bounds)
    if len(b) == 2 and not hasattr(b[0], "__len__"):
        b = [b]
    if len(b) == 1:
        b = b * n
    if len(b) != n:
        raise ValueError("Invalid input for linprog: bounds must be a single (min, max) pair "
                         "or a sequence of n (min, max) pairs")
    lb = np.empty(n)
    ub = np.empty(n)
    for j, (lo, hi) in enumerate(b):
        lb[j] = -_INF if lo is None else lo
        ub[j] = _INF if hi is None else hi
    return lb, ub


def _check_block(A, b, n, name):
    if A is None:
        if b is not None and np.size(b):
            raise ValueError(f"Invalid input for linprog: A_{name} must be given with b_{name}")
        return None, None
    A = _as_matrix(A)
    if A.shape[0] == 0:
        return None, None
    if A.shape[1] != n:
        raise ValueError(f"Invalid input for linprog: A_{name} must have exactly one column "
                         "per entry of c")
    if b is None:
        raise ValueError(f"Invalid input for linprog: b_{name} must be given with A_{name}")
    b = np.asarray(b, dtype=np.float64).ravel()
    if b.shape[0] != A.shape[0]:
        raise ValueError(f"Invalid input for linprog: b_{name} must have one entry per row "
                         f"of A_{name}")
    return A, b


def _matvec(A, x):
    return A @ x


def linprog(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None, method="highs",
            callback=None, options=None, x0=None, integrality=None):
    """Linear programming with the signature of ``scipy.optimize.linprog``.

    minimize ``c @ x`` subject to ``A_ub @ x <= b_ub``, ``A_eq @ x == b_eq`` and bounds
    (default ``0 <= x``). ``method`` is "highs", "highs-ds" or "highs-ipm"; ``callback`` and
    ``x0`` are accepted and ignored, as in SciPy's HiGHS methods.
    """
    if method not in _METHODS:
        raise ValueError(f"Unknown solver {method!r}: rusolve supports 'highs', 'highs-ds' "
                         "and 'highs-ipm'")
    c, n = _check_vec(c)
    kw, raw = _solve_options(options, _LINPROG_OPTIONS, _METHODS[method])
    A_ub, b_ub = _check_block(A_ub, b_ub, n, "ub")
    A_eq, b_eq = _check_block(A_eq, b_eq, n, "eq")
    m = Model()
    lb, ub = _linprog_bounds(bounds, n)
    integ = _integrality(integrality, n)
    m.add_vars(n, lb=lb, ub=ub, cost=c, integer=integ)
    _semi(integrality, n, m)
    n_ub = 0 if A_ub is None else A_ub.shape[0]
    n_eq = 0 if A_eq is None else A_eq.shape[0]
    if n_ub:
        m.add_matrix(A_ub, "<=", b_ub)
    if n_eq:
        m.add_matrix(A_eq, "==", b_eq)
    sol = m.solve(options=raw or None, **kw)
    status, msg = _status(sol.status)
    x = np.asarray(sol.x, dtype=np.float64) if sol.ok else None
    res = OptimizeResult(x=x, fun=sol.objective if x is not None else None)
    slack = con = None
    if x is not None:
        slack = b_ub - _matvec(A_ub, x) if n_ub else np.empty(0)
        con = b_eq - _matvec(A_eq, x) if n_eq else np.empty(0)
    duals = sol.duals if (x is not None and sol.duals is not None) else None
    rc = sol.reduced_costs if duals is not None else None
    if duals is not None:
        duals = np.asarray(duals, dtype=np.float64)
        rc = None if rc is None else np.asarray(rc, dtype=np.float64)
    res.update(
        slack=slack, con=con,
        ineqlin=OptimizeResult(residual=slack, marginals=None if duals is None else duals[:n_ub]),
        eqlin=OptimizeResult(residual=con,
                             marginals=None if duals is None else duals[n_ub:n_ub + n_eq]),
        lower=OptimizeResult(residual=None if x is None else x - _vec(lb, n),
                             marginals=None if rc is None else np.maximum(rc, 0.0)),
        upper=OptimizeResult(residual=None if x is None else _vec(ub, n) - x,
                             marginals=None if rc is None else np.minimum(rc, 0.0)),
        status=status, success=status == 0, message=f"{msg}(rusolve status: {sol.status})",
        nit=sol.iterations or 0, crossover_nit=None)
    if x is not None and integ is not False:
        res.update(mip_node_count=sol.nodes, mip_dual_bound=sol.bound, mip_gap=sol.gap)
    return res
