"""Quadratic objective (convex QP).

objective = c'x + (1/2) x'Qx. Reference values are solved by hand (KKT) or with scipy.
"""
import numpy as np
import pytest
import rusolve as rs


def _model():
    # min (1/2)(x^2 + y^2) - x - 2y, x + y <= 1, x, y >= 0  ->  x = 0, y = 1, obj = -1.5
    m = rs.Model()
    x = m.add_vars(2, lb=0, ub=10, cost=[-1.0, -2.0])
    m.add_constr(x[0] + x[1] <= 1)
    return m, x


def test_dense_q():
    m, x = _model()
    m.set_quadratic_objective([[1.0, 0.0], [0.0, 1.0]])
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(-1.5, abs=1e-6)
    assert sol.x == pytest.approx([0.0, 1.0], abs=1e-5)


def test_sparse_triplets_dict_and_numpy_agree():
    import scipy.sparse as sp

    q = np.array([[2.0, 0.5], [0.5, 1.0]])
    objs = []
    for kw in (
        {"q": q},
        {"q": sp.csr_matrix(q)},
        {"q": sp.csc_matrix(q)},
        {"q": sp.coo_matrix(q)},
        {"q": {(0, 0): 2.0, (0, 1): 0.5, (1, 1): 1.0}},  # upper triangle only: symmetrised
        {"rows": [0, 0, 1, 1], "cols": [0, 1, 0, 1], "vals": [2.0, 0.5, 0.5, 1.0]},
    ):
        m, x = _model()
        m.set_quadratic_objective(**kw)
        sol = m.solve()
        assert sol.status == "Optimal", kw
        objs.append(sol.objective)
    assert objs == pytest.approx([objs[0]] * len(objs), abs=1e-7)


def test_asymmetric_q_is_symmetrised():
    a, _ = _model()
    a.set_quadratic_objective([[1.0, 2.0], [0.0, 1.0]])
    b, _ = _model()
    b.set_quadratic_objective([[1.0, 1.0], [1.0, 1.0]])
    assert a.solve().objective == pytest.approx(b.solve().objective, abs=1e-7)


def test_clear_and_resolve():
    m, x = _model()
    m.set_quadratic_objective([[1.0, 0.0], [0.0, 1.0]])
    assert m.solve().objective == pytest.approx(-1.5, abs=1e-6)
    m.set_quadratic_objective(None)  # back to the LP: min -x - 2y -> y = 1, obj -2
    assert m.solve().objective == pytest.approx(-2.0, abs=1e-9)


def test_maximize_concave():
    # max -(1/2)(x^2 + y^2) + x + 2y, x + y <= 1  (same as the minimisation above, negated)
    m = rs.Model(sense="max")
    x = m.add_vars(2, lb=0, ub=10, cost=[1.0, 2.0])
    m.add_constr(x[0] + x[1] <= 1)
    m.set_quadratic_objective([[-1.0, 0.0], [0.0, -1.0]])
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(1.5, abs=1e-6)


def test_errors_name_the_problem():
    m, x = _model()
    with pytest.raises(ValueError, match="square"):
        m.set_quadratic_objective(np.zeros((2, 3)))
    with pytest.raises(IndexError, match="out of range"):
        m.set_quadratic_objective({(0, 5): 1.0})
    with pytest.raises(ValueError, match="finite"):
        m.set_quadratic_objective({(0, 0): float("nan")})
    with pytest.raises(ValueError, match="not both"):
        m.set_quadratic_objective([[1.0]], rows=[0], cols=[0], vals=[1.0])
    with pytest.raises(ValueError, match="go together"):
        m.set_quadratic_objective(rows=[0])


def test_nonconvex_is_refused_with_names_and_eigenvalue():
    m = rs.Model()
    x = m.add_vars(3, lb=0, ub=1, names="v")
    m.minimize([0.0, 0.0, 0.0])
    m.set_quadratic_objective([[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    with pytest.raises(ValueError, match=r"nonconvex.*eigenvalue -1\.0.*v\[?[01]"):
        m.solve()


def test_max_of_convex_is_refused_and_min_of_concave_too():
    m = rs.Model(sense="max")
    m.add_vars(2, lb=0, ub=1)
    m.set_quadratic_objective([[1.0, 0.0], [0.0, 1.0]])
    with pytest.raises(ValueError, match="maximisation"):
        m.solve()
    m.set_objective([0.0, 0.0], sense="min")  # same Q, now a minimisation: convex
    assert m.solve().status == "Optimal"


def test_miqp_is_refused():
    m = rs.Model()
    m.add_vars(2, lb=0, ub=3, integer=True)
    m.set_quadratic_objective([[1.0, 0.0], [0.0, 1.0]])
    with pytest.raises(ValueError, match="MIQP"):
        m.solve()


def test_quad_info_reports_class_and_shift():
    m, x = _model()
    m.set_quadratic_objective([[1.0, 0.0], [0.0, 1.0]])
    assert m.quad_info() is None
    m.solve()
    info = m.quad_info()
    assert info["class"] == "DCL" and info["convexity"] == "convex" and info["shift"] == 0.0
    # F'F with a diagonal perturbation at rounding level: convex after a reported shift
    rng = np.random.default_rng(1)
    f = rng.standard_normal((2, 12))
    q = f.T @ f - 1e-9 * np.eye(12)
    m = rs.Model()
    m.add_vars(12, lb=-1, ub=1, cost=rng.standard_normal(12))
    m.set_quadratic_objective(q)
    assert m.solve().status == "Optimal"
    info = m.quad_info()
    assert info["convexity"] == "convex_after_shift" and 0 < info["shift"] < 1e-7
