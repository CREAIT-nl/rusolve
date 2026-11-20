"""`rusolve.scipy.milp` / `linprog` against `scipy.optimize` (docs/design/python-api.md,
"Compatibility layers"): same status, same objective (rel 1e-6), feasible x."""
import numpy as np
import pytest
import rusolve.scipy as rss
import scipy.sparse as sp
from scipy import optimize as so
from scipy.optimize import Bounds, LinearConstraint

TOL = 1e-6


def feasible_lp(res, c, A_ub, b_ub, A_eq, b_eq, lb, ub, tol=1e-6):
    x = res.x
    assert np.all(x >= lb - tol) and np.all(x <= ub + tol)
    if A_ub is not None:
        assert np.all(A_ub @ x <= np.asarray(b_ub) + tol)
    if A_eq is not None:
        assert np.allclose(A_eq @ x, b_eq, atol=tol)
    assert np.isclose(c @ x, res.fun, rtol=1e-6, atol=1e-6)


def same_objective(a, b):
    assert a.status == b.status, (a.status, b.status, a.message, b.message)
    if a.status == 0:
        assert a.fun == pytest.approx(b.fun, rel=TOL, abs=1e-9)


# ------------------------------------------------------------------------------- result type
def test_result_is_dict_and_attrs():
    res = rss.milp([1.0, 1.0], constraints=LinearConstraint([[1, 1]], [2], [np.inf]))
    assert res["status"] == res.status == 0 and res.success is True
    assert set(res) >= {"x", "fun", "status", "success", "message", "mip_node_count",
                        "mip_dual_bound", "mip_gap"}
    assert isinstance(res.message, str)
    assert "status" in repr(res)
    with pytest.raises(AttributeError):
        res.nonexistent


def test_linprog_result_fields():
    res = rss.linprog([-1, -2], A_ub=[[1, 1], [1, 3]], b_ub=[4, 6])
    for k in ("x", "fun", "slack", "con", "status", "success", "message", "nit", "ineqlin",
              "eqlin", "lower", "upper"):
        assert k in res
    assert res.ineqlin.marginals.shape == (2,) and res.ineqlin.residual.shape == (2,)
    assert res.eqlin.marginals.shape == (0,)


# ------------------------------------------------------------------- documented milp examples
def test_milp_docs_example():
    c = -np.array([0, 1])
    A = np.array([[-1, 1], [3, 2], [2, 3]])
    b_u = np.array([1, 12, 12])
    b_l = np.full_like(b_u, -np.inf, dtype=float)
    cons = LinearConstraint(A, b_l, b_u)
    integ = np.ones_like(c)
    a = rss.milp(c=c, constraints=cons, integrality=integ)
    b = so.milp(c=c, constraints=cons, integrality=integ)
    assert list(a.x) == pytest.approx([1.0, 2.0])
    same_objective(a, b)
    a = rss.milp(c=c, constraints=cons)
    b = so.milp(c=c, constraints=cons)
    assert list(a.x) == pytest.approx([1.8, 2.8])
    same_objective(a, b)


def test_milp_tutorial_knapsack():
    values = np.array([8, 3, 12, 6, 9, 1, 7])
    sizes = np.array([3, 4, 5, 2, 6, 7, 1])
    cons = LinearConstraint(sizes[None, :], -np.inf, 12)
    kw = dict(c=-values, constraints=cons, integrality=np.ones(7), bounds=Bounds(0, 1))
    a, b = rss.milp(**kw), so.milp(**kw)
    same_objective(a, b)
    assert set(np.round(a.x)) <= {0.0, 1.0}


def test_milp_tutorial_assignment_style():
    # 3x3 assignment with equality rows, binary variables, sparse A
    rng = np.random.default_rng(5)
    cost = rng.integers(1, 20, (3, 3)).astype(float)
    rows = sp.lil_matrix((6, 9))
    for i in range(3):
        for j in range(3):
            rows[i, 3 * i + j] = 1
            rows[3 + j, 3 * i + j] = 1
    A = rows.tocsr()
    cons = LinearConstraint(A, np.ones(6), np.ones(6))
    kw = dict(c=cost.ravel(), constraints=cons, integrality=1, bounds=Bounds(0, 1))
    a, b = rss.milp(**kw), so.milp(**kw)
    same_objective(a, b)
    assert np.allclose(A @ a.x, 1.0)


# ------------------------------------------------------------------ documented linprog example
def test_linprog_docs_example():
    c = [-1, 4]
    A = [[-3, 1], [1, 2]]
    b = [6, 4]
    bounds = [(None, None), (-3, None)]
    a = rss.linprog(c, A_ub=A, b_ub=b, bounds=bounds)
    s = so.linprog(c, A_ub=A, b_ub=b, bounds=bounds)
    assert a.fun == pytest.approx(-22.0)
    assert list(a.x) == pytest.approx([10.0, -3.0])
    same_objective(a, s)
    assert a.nit >= 0


def test_linprog_marginals_match_scipy():
    c = [-1, -2]
    A = [[1, 1], [1, 3]]
    b = [4, 6]
    a = rss.linprog(c, A_ub=A, b_ub=b)
    s = so.linprog(c, A_ub=A, b_ub=b, method="highs")
    assert a.ineqlin.marginals == pytest.approx(s.ineqlin.marginals, abs=1e-8)
    assert a.slack == pytest.approx(s.slack, abs=1e-8)
    assert a.lower.marginals == pytest.approx(s.lower.marginals, abs=1e-8)
    assert a.upper.marginals == pytest.approx(s.upper.marginals, abs=1e-8)


def test_linprog_eq_marginals_and_con():
    a = rss.linprog([1, 2], A_eq=[[1, 1]], b_eq=[3], bounds=(0, 2))
    s = so.linprog([1, 2], A_eq=[[1, 1]], b_eq=[3], bounds=(0, 2))
    same_objective(a, s)
    assert a.eqlin.marginals == pytest.approx(s.eqlin.marginals, abs=1e-8)
    assert a.con == pytest.approx([0.0], abs=1e-8)


# --------------------------------------------------------------------------- bounds forms
BOUND_CASES = [
    None,                              # default (0, None)
    (0, 4),                            # single pair
    (None, 4),
    [(0, 1), (-2, 3), (None, None)],   # per variable pairs
    [(0, None)] * 3,
]


@pytest.mark.parametrize("bounds", BOUND_CASES)
def test_linprog_bounds_forms(bounds):
    c = np.array([1.0, -1.0, 0.5])
    A = np.array([[1.0, 1.0, 1.0]])
    kw = dict(A_ub=A, b_ub=[5.0], A_eq=[[1.0, 0.0, -1.0]], b_eq=[0.0], bounds=bounds)
    a, s = rss.linprog(c, **kw), so.linprog(c, **kw)
    same_objective(a, s)


def test_linprog_bounds_object():
    # scipy.optimize.linprog itself rejects Bounds objects: compare with the pair-list form
    a = rss.linprog([-1.0, -1.0], A_ub=[[1.0, 1.0]], b_ub=[3.0], bounds=Bounds([0, 1], [2, 5]))
    s = so.linprog([-1.0, -1.0], A_ub=[[1.0, 1.0]], b_ub=[3.0], bounds=[(0, 2), (1, 5)])
    same_objective(a, s)


@pytest.mark.parametrize("bounds", [None, Bounds(0, 1), Bounds([0, -1], [2, 3]), (0, 3)])
def test_milp_bounds_forms(bounds):
    kw = dict(c=[-1.0, -2.0], constraints=LinearConstraint([[1, 1]], -np.inf, 2.5),
              integrality=[1, 1], bounds=bounds)
    a, s = rss.milp(**kw), so.milp(**kw)
    same_objective(a, s)


# ------------------------------------------------------------------- constraints forms
def test_milp_constraint_forms():
    c = [-1.0, -1.0, -1.0]
    A1 = np.array([[1.0, 1.0, 0.0]])
    A2 = sp.csr_matrix([[0.0, 1.0, 1.0], [1.0, 0.0, 1.0]])
    c1 = LinearConstraint(A1, -np.inf, 3.0)
    c2 = LinearConstraint(A2, [1.0, -np.inf], [np.inf, 4.0])
    ub = Bounds(0, 10)
    cases = [c1, [c1], [c1, c2], (c1, c2)]
    for cons in cases:
        a = rss.milp(c, constraints=cons, bounds=ub, integrality=1)
        s = so.milp(c, constraints=cons, bounds=ub, integrality=1)
        same_objective(a, s)
    # (A, lb, ub) tuples
    a = rss.milp(c, constraints=[(A1, -np.inf, 3.0), (A2, [1.0, -np.inf], [np.inf, 4.0])],
                 bounds=ub, integrality=1)
    same_objective(a, so.milp(c, constraints=[c1, c2], bounds=ub, integrality=1))


def test_milp_no_constraints_and_ranged_rows():
    # bounds only
    a = rss.milp([1.0, -1.0], bounds=Bounds(0, 3))
    same_objective(a, so.milp([1.0, -1.0], bounds=Bounds(0, 3)))
    # ranged, equality and one-sided rows mixed in one block
    A = np.array([[1.0, 1.0], [1.0, -1.0], [2.0, 1.0], [1.0, 0.0]])
    lo = np.array([1.0, 0.0, -np.inf, -np.inf])
    hi = np.array([4.0, 0.0, 6.0, np.inf])  # last row is free
    cons = LinearConstraint(A, lo, hi)
    for c in ([-1.0, -1.0], [1.0, 1.0], [1.0, -3.0]):
        a = rss.milp(c, constraints=cons)
        same_objective(a, so.milp(c, constraints=cons))


# ---------------------------------------------------------------------- random problems
def random_lp(seed, m=8, n=6, sparse=False):
    rng = np.random.default_rng(seed)
    A = rng.integers(-3, 6, (m, n)).astype(float)
    A[rng.random((m, n)) < 0.4] = 0.0
    x0 = rng.uniform(0, 3, n)
    b = A @ x0 + rng.uniform(0, 2, m)
    c = rng.integers(-5, 6, n).astype(float)
    ub = np.full(n, 8.0)
    return c, (sp.csr_matrix(A) if sparse else A), b, ub


@pytest.mark.parametrize("seed", range(12))
@pytest.mark.parametrize("sparse", [False, True])
def test_random_lp(seed, sparse):
    c, A, b, ub = random_lp(seed, sparse=sparse)
    bounds = list(zip(np.zeros_like(ub), ub))
    a = rss.linprog(c, A_ub=A, b_ub=b, bounds=bounds)
    s = so.linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
    same_objective(a, s)
    if a.status == 0:
        feasible_lp(a, c, A, b, None, None, 0.0, ub)


@pytest.mark.parametrize("seed", range(8))
def test_random_lp_equality_only(seed):
    rng = np.random.default_rng(100 + seed)
    m, n = 4, 7
    A = rng.integers(-2, 4, (m, n)).astype(float)
    x0 = rng.uniform(0.5, 2, n)
    b = A @ x0
    c = rng.uniform(0, 3, n)
    a = rss.linprog(c, A_eq=A, b_eq=b, bounds=(0, 6))
    s = so.linprog(c, A_eq=A, b_eq=b, bounds=(0, 6), method="highs")
    same_objective(a, s)
    if a.status == 0:
        feasible_lp(a, c, None, None, A, b, 0.0, 6.0)


@pytest.mark.parametrize("seed", range(6))
def test_random_lp_free_variables(seed):
    rng = np.random.default_rng(200 + seed)
    m, n = 5, 4
    A = rng.integers(-3, 4, (m, n)).astype(float)
    b = A @ rng.uniform(-1, 1, n) + rng.uniform(0.1, 1, m)
    c = rng.uniform(-1, 1, n)
    # box the problem through extra rows so it stays bounded
    A2 = np.vstack([A, np.eye(n), -np.eye(n)])
    b2 = np.concatenate([b, np.full(2 * n, 5.0)])
    bounds = [(None, None)] * n
    a = rss.linprog(c, A_ub=A2, b_ub=b2, bounds=bounds)
    s = so.linprog(c, A_ub=A2, b_ub=b2, bounds=bounds, method="highs")
    same_objective(a, s)


@pytest.mark.parametrize("seed", range(12))
@pytest.mark.parametrize("sparse", [False, True])
def test_random_mip(seed, sparse):
    rng = np.random.default_rng(300 + seed)
    m, n = 5, 8
    A = rng.integers(0, 6, (m, n)).astype(float)
    b = A.sum(axis=1) * rng.uniform(0.3, 0.6, m)
    c = -rng.integers(1, 10, n).astype(float)
    integrality = (rng.random(n) < 0.7).astype(int)
    M = sp.csr_matrix(A) if sparse else A
    cons = LinearConstraint(M, -np.inf, b)
    kw = dict(c=c, constraints=cons, integrality=integrality, bounds=Bounds(0, 4))
    a, s = rss.milp(**kw), so.milp(**kw)
    same_objective(a, s)
    if a.status == 0:
        assert np.all(A @ a.x <= b + 1e-6)
        assert np.all(a.x >= -1e-6) and np.all(a.x <= 4 + 1e-6)
        ints = integrality.astype(bool)
        assert np.allclose(a.x[ints], np.round(a.x[ints]), atol=1e-6)
        assert a.mip_gap is None or a.mip_gap <= 1e-4 + 1e-9


@pytest.mark.parametrize("seed", range(4))
def test_random_mip_equality_rows(seed):
    rng = np.random.default_rng(400 + seed)
    n = 7
    A = rng.integers(0, 4, (2, n)).astype(float)
    xs = rng.integers(0, 3, n)
    b = A @ xs
    c = rng.integers(1, 9, n).astype(float)
    cons = LinearConstraint(A, b, b)
    kw = dict(c=c, constraints=cons, integrality=1, bounds=Bounds(0, 5))
    same_objective(rss.milp(**kw), so.milp(**kw))


def test_linprog_with_integrality():
    c, A, b, ub = random_lp(7)
    kw = dict(A_ub=A, b_ub=b, bounds=list(zip(np.zeros_like(ub), ub)), integrality=1)
    a, s = rss.linprog(c, **kw), so.linprog(c, **kw)
    same_objective(a, s)
    if a.status == 0:
        assert np.allclose(a.x, np.round(a.x), atol=1e-6)
        assert "mip_gap" in a


# ----------------------------------------------------------------------- status handling
def test_infeasible_lp():
    kw = dict(A_ub=[[1.0]], b_ub=[1.0], A_eq=[[1.0]], b_eq=[3.0])
    a, s = rss.linprog([1.0], **kw), so.linprog([1.0], **kw)
    assert a.status == s.status == 2
    assert a.success is False and a.x is None and a.fun is None


def test_infeasible_milp():
    cons = LinearConstraint([[2.0]], [1.0], [1.0])
    a = rss.milp([1.0], constraints=cons, integrality=1)
    s = so.milp([1.0], constraints=cons, integrality=1)
    assert a.status == s.status == 2
    assert a.x is None and a.fun is None and a.success is False


def test_unbounded_lp():
    kw = dict(A_ub=[[1.0, -1.0]], b_ub=[1.0], bounds=(0, None))
    a, s = rss.linprog([-1.0, -1.0], **kw), so.linprog([-1.0, -1.0], **kw)
    assert a.status == 3
    assert s.status in (3, 4)
    assert a.success is False


def test_unbounded_milp():
    kw = dict(c=[-1.0, -1.0], constraints=LinearConstraint([[1.0, -1.0]], -np.inf, 1.0),
              integrality=1)
    a, s = rss.milp(**kw), so.milp(**kw)
    assert a.status in (3, 4) and s.status in (3, 4)


def test_unbounded_variables_with_free_bounds():
    # min x s.t. x >= -3 (a row), x free
    a = rss.milp([1.0], bounds=Bounds(-np.inf, np.inf),
                 constraints=LinearConstraint([[1.0]], -3.0, np.inf))
    assert a.x[0] == pytest.approx(-3.0)
    a = rss.linprog([1.0], A_ub=[[-1.0]], b_ub=[3.0], bounds=(None, None))
    assert a.x[0] == pytest.approx(-3.0)


# ----------------------------------------------------------------------------- options
def test_milp_options():
    kw = dict(c=[-1.0, -1.0], constraints=LinearConstraint([[1, 1]], -np.inf, 3.5),
              integrality=1)
    a = rss.milp(**kw, options={"disp": False, "presolve": True, "time_limit": 30,
                                "mip_rel_gap": 1e-6, "node_limit": 1000})
    assert a.status == 0 and a.fun == pytest.approx(-3.0)
    # SciPy warns and passes unknown options to HiGHS; rusolve's solver rejects them
    with pytest.warns(RuntimeWarning), pytest.raises(ValueError):
        rss.milp(**kw, options={"no_such_option_xyz": 1})


def test_options_dict_not_mutated():
    opts = {"disp": False, "time_limit": 5}
    rss.milp([1.0], options=opts)
    rss.linprog([1.0], options=opts)
    assert opts == {"disp": False, "time_limit": 5}


def test_linprog_options_and_methods():
    kw = dict(A_ub=[[1.0, 1.0]], b_ub=[2.0])
    ref = so.linprog([-1.0, -2.0], **kw).fun
    for method in ("highs", "highs-ds", "highs-ipm"):
        a = rss.linprog([-1.0, -2.0], **kw, method=method,
                        options={"presolve": False, "time_limit": 20, "disp": False})
        assert a.status == 0 and a.fun == pytest.approx(ref, rel=1e-6)
    with pytest.raises(ValueError):
        rss.linprog([1.0], method="interior-point")
    # callback and x0 are accepted
    a = rss.linprog([-1.0, -2.0], **kw, callback=lambda r: None, x0=[0.0, 0.0])
    assert a.status == 0


# ----------------------------------------------------------------------------- validation
def test_bad_inputs():
    with pytest.raises(ValueError):
        rss.milp([])
    with pytest.raises(ValueError):
        rss.milp([1.0, 2.0], constraints=LinearConstraint([[1.0, 1.0, 1.0]], 0, 1))
    with pytest.raises(ValueError):
        rss.milp([1.0], integrality=[7])
    with pytest.raises(ValueError):
        rss.linprog([1.0, 1.0], A_ub=[[1.0]], b_ub=[1.0])
    with pytest.raises(ValueError):
        rss.linprog([1.0], A_ub=[[1.0]])  # b_ub missing
    with pytest.raises(ValueError):
        rss.linprog([1.0, 1.0], bounds=[(0, 1)] * 3)


def test_scipy_is_not_imported_by_the_shim():
    import subprocess
    import sys
    code = ("import sys; import rusolve.scipy as r; "
            "r.milp([1.0]); "
            "sys.exit(1 if any(m == 'scipy' or m.startswith('scipy.') for m in sys.modules) else 0)")
    # numpy is imported; scipy must not be (A is duck-typed)
    p = subprocess.run([sys.executable, "-c", code], capture_output=True)
    assert p.returncode == 0, p.stderr.decode()


def test_sparse_never_densified():
    # 200k x 200k sparse identity: a dense copy would need 320 GB
    n = 200_000
    A = sp.identity(n, format="csr")
    res = rss.milp(np.ones(n), constraints=LinearConstraint(A, 1.0, np.inf))
    assert res.status == 0 and res.fun == pytest.approx(float(n))
