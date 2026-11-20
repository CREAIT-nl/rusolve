"""Core tests of the native `rusolve.Model` API (docs/design/python-api.md: "The API",
"Expression layer", "Semantics")."""

import array
import itertools
import threading
import time

import numpy as np
import pytest

import rusolve as rs


def lp3():
    """min x0 + x1 + x2  s.t.  x0 + 2 x1 + x2 >= 3,  0 <= x <= 4."""
    m = rs.Model()
    x = m.add_vars(3, ub=4, cost=1.0)
    m.add_constr(x[0] + 2 * x[1] + x[2] >= 3)
    return m, x


# ---------------------------------------------------------------- exports and Vars


def test_exports():
    for name in ["Model", "Var", "Vars", "LinExpr", "Constraint", "Solution", "quicksum"]:
        assert hasattr(rs, name)
    assert rs.inf == float("inf")


def test_vars_block():
    m = rs.Model()
    x = m.add_vars(5)
    assert len(x) == 5 and x.start == 0 and x.stop == 5
    assert x.indices == range(0, 5)
    assert x[-1].index == 4 and int(x[2]) == 2
    assert [v.index for v in x] == [0, 1, 2, 3, 4]
    assert [v.index for v in x[1:4]] == [1, 2, 3]
    assert [v.index for v in x[::2]] == [0, 2, 4]
    assert len(x[3:1]) == 0
    with pytest.raises(IndexError):
        x[5]
    y = m.add_vars(2)
    assert y.start == 5 and m.num_vars == 7
    assert repr(x[1]) == "Var(1)"
    assert list(range(10))[x[3]] == 3  # __index__


def test_var_hash_and_equality():
    m = rs.Model()
    x = m.add_vars(3)
    d = {x[0]: "a", x[1]: "b"}
    assert d[x[0]] == "a" and d[m.var(1)] == "b"
    assert x[2] not in d
    assert x[0] in [x[1], x[0]]
    assert x[0] == x[0] and not (x[0] == x[1])
    assert x[0] != x[1] and not (x[0] != x[0])
    other = rs.Model().add_vars(1)
    assert not (x[0] == other[0]) and {x[0]: 1}.get(other[0]) is None
    c = x[0] <= x[1]
    with pytest.raises(TypeError, match="between"):
        bool(c)
    with pytest.raises(TypeError, match="between"):
        if x[0] + 1 == 3:
            pass
    with pytest.raises(TypeError, match="strict"):
        x[0] < 1


def test_linexpr_arithmetic():
    m = rs.Model()
    x = m.add_vars(3)
    e = 2 * x[0] - x[1] / 2 + 3
    assert isinstance(e, rs.LinExpr)
    assert e.constant == 3
    assert sorted(e.terms()) == [(0, 2.0), (1, -0.5)]
    f = -e
    assert f.constant == -3 and sorted(f.terms()) == [(0, -2.0), (1, 0.5)]
    g = 1 - x[2]
    assert g.constant == 1 and g.terms() == [(2, -1.0)]
    h = rs.LinExpr()
    h += x[0]
    h += 2 * x[1]
    h -= 1
    assert h.constant == -1 and sorted(h.terms()) == [(0, 1.0), (1, 2.0)]
    h += h
    assert h.constant == -2 and len(h) == 4
    # e is not modified by expressions built from it
    k = e + x[2]
    assert len(e) == 2 and len(k) == 3


def test_sum_quicksum_dot():
    m = rs.Model()
    x = m.add_vars(4)
    s = sum(x)
    assert sorted(s.terms()) == [(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0)]
    assert sorted(x.sum().terms()) == sorted(s.terms())
    q = rs.quicksum([x[0], 2 * x[1], 5, (3.0, x[2]), x[3:]])
    assert q.constant == 5
    assert sorted(q.terms()) == [(0, 1.0), (1, 2.0), (2, 3.0), (3, 1.0)]
    q2 = rs.quicksum(c * v for c, v in zip([1, 2, 3, 4], x))
    assert sorted(q2.terms()) == [(0, 1.0), (1, 2.0), (2, 3.0), (3, 4.0)]
    d = x.dot(np.array([1.0, 0.0, 2.0, 3.0]))
    assert sorted(d.terms()) == [(0, 1.0), (1, 0.0), (2, 2.0), (3, 3.0)] or sorted(
        d.terms()
    ) == [(0, 1.0), (2, 2.0), (3, 3.0)]
    with pytest.raises(ValueError):
        x.dot([1.0, 2.0])


def test_constraint_objects():
    m = rs.Model()
    x = m.add_vars(2)
    c = x[0] + 2 * x[1] + 1 <= 3
    assert (c.lo, c.hi) == (-rs.inf, 2.0)
    c = 3 >= x[0]
    assert (c.lo, c.hi) == (-rs.inf, 3.0)
    c = x[0] + 1 == x[1] - 2
    assert (c.lo, c.hi) == (-3.0, -3.0)
    c = (x[0] + x[1]).between(1, 4)
    assert (c.lo, c.hi) == (1.0, 4.0)
    c = x[0].between(-1, 1)
    assert (c.lo, c.hi) == (-1.0, 1.0)


# ---------------------------------------------------------------- rows


def test_every_row_form():
    m = rs.Model()
    x = m.add_vars(3, ub=10)
    rows = [
        m.add_constr({x[0]: 1.0, x[1]: 2.0}, "<=", 4),
        m.add_constr({0: 1.0, 2: 1.0}, ">=", 1),
        m.add_constr(([0, 1], [1.0, 1.0]), "==", 2),
        m.add_constr((np.array([1, 2], dtype=np.int32), np.array([1.0, -1.0])), "L", 3),
        m.add_constr((np.array([0, 2], dtype=np.int64), [2.0, 2.0]), "G", 0),
        m.add_constr((array.array("i", [0]), array.array("d", [1.0])), "E", 1),
        m.add_constr((array.array("q", [1]), array.array("d", [1.0])), "<=", 9),
        m.add_constr(x[0] + x[1] - 1, "<=", 5),  # constant moves to the bound
        m.add_constr(x[2], ">=", 0.5),
        m.add_constr(x[0] + x[1] <= 7),
        m.add_constr(x[1] >= x[2]),
        m.add_constr(x[0] + x[1] + x[2], lo=1, hi=8),
        m.add_constr({x[0]: 1}, lo=0),
        m.add_constr({x[1]: 1}, hi=10),
        m.add_constr((x[0] - x[2]).between(-5, 5)),
        m.add_constr({x[0]: 1.0, x[1]: 1.0}, "=<", 100),
    ]
    assert rows == list(range(16))
    assert m.num_constrs == 16
    sol = m.solve()
    assert sol.status == "Optimal"
    act = sol.row_activity
    x0, x1, x2 = sol.x
    assert act[2] == pytest.approx(2.0) and x0 + x1 == pytest.approx(2.0)
    assert x0 == pytest.approx(1.0)  # row 5: x0 == 1
    assert act[7] == pytest.approx(x0 + x1)
    assert 1 - 1e-9 <= act[11] <= 8 + 1e-9


def test_dict_subclass_and_row_iterables():
    from collections import OrderedDict

    m = rs.Model()
    x = m.add_vars(3, ub=1)
    rr = m.add_constrs([OrderedDict([(x[0], 1.0), (x[1], 1.0)]), {2: 1.0}], "<=", 1)
    assert rr == range(0, 2)
    rr = m.add_constrs((x[k] <= 0.5 for k in range(3)))
    assert rr == range(2, 5)
    rr = m.add_constrs(({k: 1.0} for k in range(3)), [">=", ">=", "<="], np.zeros(3))
    assert rr == range(5, 8)
    rr = m.add_constrs([{0: 1.0}, {1: 1.0}], lo=[0, 0], hi=np.array([1.0, 2.0]))
    assert rr == range(8, 10)
    assert m.add_constrs([], "<=", 1) == range(10, 10)
    assert m.solve().status == "Optimal"


def test_dict_rows_fast_and_slow_entries():
    """The raw PyDict_Next path (exact int keys, exact float/int values) and its fallbacks
    (Var keys, numpy scalars, bool, out-of-range ints) give the same rows and errors."""
    m = rs.Model()
    x = m.add_vars(3, ub=10)
    m.add_constr({0: 2, 1: 1.5, x[2]: np.float64(1.0)}, "<=", 8)  # int, float, numpy values
    m.add_constr({np.int64(0): 1.0, True: 1}, ">=", 1)  # numpy and bool keys
    m.maximize(x[0] + x[1] + x[2])
    sol = m.solve()
    assert sol.status == "Optimal"
    assert m.num_nonzeros == 5
    with pytest.raises(TypeError, match="variable"):
        m.add_constr({2**70: 1.0}, "<=", 1)
    with pytest.raises(OverflowError):
        m.add_constr({0: 10**400}, "<=", 1)
    with pytest.raises(TypeError, match="number"):
        m.add_constr({0: "1"}, "<=", 1)
    assert m.num_constrs == 2


def test_dict_row_mutated_by_value_conversion():
    """A value whose __float__ mutates the dict being read must not crash the reader."""
    m = rs.Model()
    m.add_vars(64, ub=1)
    row = {}

    class Evil:
        def __float__(self):
            row.clear()
            for k in range(32, 64):
                row[k] = 1.0
            return 1.0

    row.update({k: 1.0 for k in range(32)})
    row[5] = Evil()
    try:
        m.add_constr(row, "<=", 1)
    except (RuntimeError, ValueError, TypeError):
        pass
    m.solve()


def test_duplicates_and_zeros_merged():
    m = rs.Model()
    x = m.add_vars(3, ub=10)
    m.add_constr(([0, 0, 1, 2], [1.0, 2.0, 0.0, 1.0]), "<=", 6)  # 3 x0 + x2
    m.add_constr(x[0] + x[0] + x[1] - x[1], ">=", 1)  # 2 x0 (x1 cancels)
    assert m.num_nonzeros == 3
    m.maximize(x[0])
    sol = m.solve()
    assert sol.objective == pytest.approx(2.0)


def test_numpy_list_array_inputs_for_vars():
    for lb, ub, cost in [
        ([0, 1, 2], [5, 5, 5], [1, 1, 1]),
        (np.array([0.0, 1, 2]), np.array([5.0, 5, 5]), np.ones(3)),
        (array.array("d", [0, 1, 2]), array.array("d", [5, 5, 5]), array.array("d", [1, 1, 1])),
        ((0, 1, 2), 5, 1.0),
        (np.array([0, 1, 2], dtype=np.int64), np.float32(5), np.ones(3, dtype=np.float32)),
        (np.arange(6.0)[::2] / 2, 5, 1),  # non-contiguous buffer
    ]:
        m = rs.Model()
        x = m.add_vars(3, lb=lb, ub=ub, cost=cost)
        sol = m.solve()
        assert sol.status == "Optimal"
        assert sol.objective == pytest.approx(3.0)
        assert list(sol[x]) == pytest.approx([0, 1, 2])


def test_add_var_and_integer_flags():
    m = rs.Model()
    a = m.add_var(lb=-1, ub=2.5, cost=-1, integer=True)
    b = m.add_var(binary=True, cost=-1)
    c = m.add_var(lb=0, ub=1, binary=True, cost=-1)  # explicit [0, 1] is allowed
    xs = m.add_vars(2, ub=1.5, integer=[True, False], cost=-1)
    sol = m.solve()
    assert sol[a] == pytest.approx(2) and sol[b] == pytest.approx(1) and sol[c] == pytest.approx(1)
    assert list(sol[xs]) == pytest.approx([1, 1.5])



def test_binary_default_bounds_explicit_bounds_apply():
    m = rs.Model(sense="max")
    a = m.add_var(binary=True, cost=1)
    b = m.add_var(binary=True, ub=3, cost=1)  # explicit ub still applies: integer in [0, 3]
    c = m.add_var(binary=True, lb=1, cost=-1)  # fixed at 1
    xs = m.add_vars(2, ub=1, binary=True, cost=1)
    sol = m.solve()
    assert (sol[a], sol[b], sol[c]) == pytest.approx((1, 3, 1))
    assert list(sol[xs]) == pytest.approx([1, 1])


def test_add_matrix_forms():
    sp = pytest.importorskip("scipy.sparse")
    dense = np.array([[1.0, 2.0, 0.0], [0.0, 1.0, 1.0]])
    forms = [
        sp.csr_matrix(dense),
        sp.csc_matrix(dense),
        sp.coo_matrix(dense),
        sp.lil_matrix(dense),
        sp.csr_array(dense),
        dense,
        dense.tolist(),
        np.asfortranarray(dense),
    ]
    for A in forms:
        m = rs.Model()
        x = m.add_vars(3, ub=10, cost=-1)
        rr = m.add_matrix(A, "<=", np.array([4.0, 3.0]), names="m")
        assert rr == range(0, 2), type(A)
        assert m.constr_name(1) == "m[1]"
        assert m.num_nonzeros == 4, type(A)
        sol = m.solve()
        assert sol.objective == pytest.approx(-7.0), type(A)  # x = (4, 0, 3)
    # vars maps the matrix columns onto chosen variables; lo/hi ranged rows
    m = rs.Model()
    x = m.add_vars(4, ub=10)
    m.add_matrix(sp.csr_matrix([[1.0, 1.0]]), vars=x[2:], lo=2, hi=3)
    m.add_matrix(np.array([[1.0, 1.0]]), "==", 1, vars=[x[0], x[3]])
    m.minimize(x.sum())
    sol = m.solve()
    assert sol.objective == pytest.approx(2.0)
    assert sol.x[0] + sol.x[3] == pytest.approx(1.0)
    with pytest.raises(ValueError, match="columns"):
        m.add_matrix(dense, "<=", 1, vars=x[:2])


# ---------------------------------------------------------------- errors


def test_error_messages():
    m = rs.Model()
    x = m.add_vars(2)
    other = rs.Model().add_vars(1)
    with pytest.raises(ValueError, match="NaN"):
        m.add_constr({x[0]: float("nan")}, "<=", 1)
    with pytest.raises(ValueError, match="(?i)infinite|inf"):
        m.add_constr({x[0]: float("inf")}, "<=", 1)
    with pytest.raises(ValueError, match='"cap"'):
        m.add_constr({x[1]: float("nan")}, "<=", 1, name="cap")
    with pytest.raises(IndexError, match="5"):
        m.add_constr({5: 1.0}, "<=", 1)
    with pytest.raises(IndexError):
        m.add_constr(([0, 7], [1.0, 1.0]), "<=", 1)
    with pytest.raises(ValueError, match="(?i)model"):
        m.add_constr(x[0] + other[0] <= 1)
    with pytest.raises(ValueError, match="(?i)model"):
        m.add_constr({other[0]: 1.0}, "<=", 1)
    with pytest.raises(TypeError):
        m.add_constr("x0 <= 1")
    with pytest.raises(ValueError, match="sense"):
        m.add_constr({0: 1.0}, "<", 1)
    with pytest.raises(ValueError, match="bounds"):
        m.add_constr({0: 1.0})
    with pytest.raises(ValueError):
        m.add_constr(x[0] <= 1, "<=", 2)
    with pytest.raises(ValueError, match="lb"):
        m.add_vars(3, lb=[0.0, 1.0])
    with pytest.raises(ValueError, match="names"):
        m.add_vars(2, names=["a"])
    with pytest.raises(ValueError, match="(?i)option"):
        m.solve(options={"no_such_option": 1})
    with pytest.raises(ValueError, match="sense"):
        rs.Model(sense="sideways")
    assert m.num_constrs == 0  # nothing failed half-way
    assert m.num_vars == 2


def test_failed_set_calls_change_nothing():
    m = rs.Model()
    x = m.add_vars(3, ub=5, cost=1)
    m.add_constr(x.sum() >= 3)
    assert m.solve().objective == pytest.approx(3)
    y = m.add_vars(2, ub=5, cost=1)  # unsynced columns
    with pytest.raises(ValueError, match="NaN"):
        m.set_bounds([x[0], y[0], y[1]], [1, 1, float("nan")], 5)
    with pytest.raises(ValueError, match="finite"):
        m.set_cost([y[0], y[1]], [2.0, float("inf")])
    with pytest.raises(ValueError, match="NaN"):
        m.set_row_bounds([0], float("nan"), 1)
    with pytest.raises(ValueError, match="(?i)model"):
        m.set_bounds(rs.Model().add_vars(1), 0, 1)
    with pytest.raises(ValueError, match="NaN"):
        m.add_vars(2, lb=[0, float("nan")])
    assert m.num_vars == 5
    sol = m.solve()
    assert sol.objective == pytest.approx(3)  # neither the bounds nor the costs changed
    assert list(sol[y]) == pytest.approx([0, 0])


def test_add_constrs_is_all_or_nothing():
    m = rs.Model()
    x = m.add_vars(2)
    m.add_constr(x[0] <= 1)
    with pytest.raises(IndexError):
        m.add_constrs([{0: 1.0}, {1: 1.0}, {9: 1.0}], "<=", 1)
    assert m.num_constrs == 1 and m.num_nonzeros == 1
    with pytest.raises(ValueError, match="rhs"):
        m.add_constrs([{0: 1.0}, {1: 1.0}], "<=", [1.0])
    assert m.num_constrs == 1
    assert m.solve().status == "Optimal"


# ---------------------------------------------------------------- objective and status


def test_objective_forms():
    m = rs.Model()
    x = m.add_vars(3, lb=1, ub=2)
    for obj, expect in [
        ({x[0]: 1.0, 2: 2.0}, 3.0),
        (x[0] + 2 * x[1] + 10, 13.0),
        (x[2], 1.0),
        (([0, 1], [1.0, 1.0]), 2.0),
        (np.array([1.0, 1.0, 1.0]), 3.0),
        ([0, 0, 3], 3.0),
    ]:
        m.set_objective(obj)
        assert m.solve().objective == pytest.approx(expect)
    m.maximize(x.sum())
    assert m.sense == "max"
    assert m.solve().objective == pytest.approx(6.0)
    m.set_objective({0: 1.0}, sense="min")
    assert m.solve().objective == pytest.approx(1.0)
    with pytest.raises(ValueError, match="3"):
        m.set_objective([1.0, 2.0])


def test_maximize_lp_duals():
    # max 3x + 2y  s.t.  x + y <= 4,  x + 3y <= 6: optimum (4, 0), only row 0 binds
    m = rs.Model(sense="max")
    x, y = m.add_var(cost=3), m.add_var(cost=2)
    r0 = m.add_constr(x + y <= 4)
    r1 = m.add_constr(x + 3 * y <= 6)
    sol = m.solve()
    assert sol.status == "Optimal" and sol.ok
    assert sol.objective == pytest.approx(12.0)
    assert sol.bound == pytest.approx(12.0) and sol.gap == 0.0
    assert sol[x] == pytest.approx(4) and sol[y] == pytest.approx(0)
    assert sol.duals is not None and len(sol.duals) == 2
    assert abs(sol.duals[r0]) == pytest.approx(3.0) and sol.duals[r1] == pytest.approx(0.0)
    assert sol.reduced_costs is not None and len(sol.reduced_costs) == 2
    assert sol.nodes == 0 and sol.iterations >= 1 and sol.time >= 0
    assert isinstance(sol.x, np.ndarray) and sol.x.dtype == np.float64
    assert "Optimal" in repr(sol)


def test_infeasible_and_unbounded():
    m = rs.Model()
    z = m.add_var()
    m.add_constr(z >= 2)
    m.add_constr(z <= 1)
    sol = m.solve()
    assert sol.status == "Infeasible"
    assert not sol.ok and sol.x is None and sol.objective is None and sol.duals is None
    with pytest.raises(ValueError, match="Infeasible"):
        sol[z]

    m = rs.Model()
    z = m.add_var(lb=-rs.inf)
    m.minimize(z)
    sol = m.solve()
    assert sol.status in ("Unbounded", "InfeasibleOrUnbounded")
    assert m.solve(presolve=False).status == "Unbounded"


def test_mip_knapsack():
    values = [10, 13, 7, 8, 5]
    weights = [5, 6, 4, 3, 2]
    best = max(
        sum(v for v, s in zip(values, pick) if s)
        for pick in itertools.product([0, 1], repeat=5)
        if sum(w for w, s in zip(weights, pick) if s) <= 10
    )
    assert best == 23
    m = rs.Model()
    x = m.add_vars(5, binary=True, names="take")
    m.maximize(x.dot(values))
    m.add_constr(x.dot(weights) <= 10, name="capacity")
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(23)
    assert sol.bound == pytest.approx(23) and sol.gap == pytest.approx(0, abs=1e-9)
    assert sol.duals is None and sol.reduced_costs is None
    items = {f"item{k}": x[k] for k in range(5)}
    assert sorted(sol.selected(items)) == ["item0", "item3", "item4"]
    assert sol.values(items, above=0.5) == pytest.approx({"item0": 1, "item3": 1, "item4": 1})
    assert sol.values([x[0], x[1]]) == pytest.approx([1, 0])
    assert sol[x.dot(weights)] == pytest.approx(10)
    assert sol[x.dot(weights) + 5] == pytest.approx(15)  # the constant is included
    assert sol[{0: 1.0, 3: 2.0}] == pytest.approx(3)
    assert sol[0] == pytest.approx(1)


# ---------------------------------------------------------------- persistence


def fresh_after_change():
    m = rs.Model(sense="max")
    x = m.add_vars(3, ub=[4, 3, 5], cost=[3.0, 2.0, 4.0])
    m.add_constr(x[0] + x[1] + x[2] <= 6)
    m.add_constr(2 * x[0] + x[1] + 3 * x[2] <= 15)
    return m, x


def test_resolve_after_changes_is_warm():
    # fresh solves without presolve, so both count simplex iterations from a basis
    m, x = fresh_after_change()
    s0 = m.solve()
    assert s0.status == "Optimal"
    x0 = s0.x.copy()

    m.set_bounds(x[2], 0, 2)
    s1 = m.solve()
    f, fx = fresh_after_change()
    f.set_bounds(fx[2], 0, 2)
    f1 = f.solve(presolve=False)
    assert s1.objective == pytest.approx(f1.objective)
    assert s1.iterations <= f1.iterations

    m.set_cost(x[0], 10.0)
    s2 = m.solve()
    f.set_cost(fx[0], 10.0)
    assert s2.objective == pytest.approx(f.solve().objective)

    m.add_constr(x[0] <= 1)
    s3 = m.solve()
    g, gx = fresh_after_change()
    g.set_bounds(gx[2], 0, 2)
    g.set_cost(gx[0], 10.0)
    g.add_constr(gx[0] <= 1)
    g3 = g.solve(presolve=False)
    assert s3.objective == pytest.approx(g3.objective)
    assert s3.iterations <= g3.iterations
    # snapshot semantics: earlier solutions do not change
    assert np.array_equal(s0.x, x0)


def test_changes_on_synced_and_new_entities():
    m = rs.Model()
    x = m.add_vars(2, ub=10, cost=1)
    r = m.add_constr(x[0] + x[1] >= 4)
    assert m.solve().objective == pytest.approx(4)
    y = m.add_var(ub=10, cost=0.5)  # new column, then a coefficient in a synced row
    m.set_coef(r, y, 1.0)
    assert m.solve().objective == pytest.approx(2)
    m.set_coef(r, y, 0.0)
    m.set_row_bounds(r, 6, rs.inf)
    assert m.solve().objective == pytest.approx(6)
    r2 = m.add_constr({x[0]: 1.0}, ">=", 1)
    m.set_row_bounds(r2, 2, rs.inf)  # unsynced row
    m.set_coef(r2, x[1], 1.0)  # unsynced row, synced column
    s = m.solve()
    assert s.objective == pytest.approx(6) and s.row_activity[r2] >= 2 - 1e-9
    m.set_bounds([x[0], x[1]], [0.5, 0], [0.5, 10])
    m.set_integer(x[1])
    s = m.solve()
    assert s.objective == pytest.approx(6.5) and s[x[1]] == pytest.approx(6)
    m.set_integer(x[1], False)
    m.set_bounds(x, 0, 10)
    assert m.solve().objective == pytest.approx(6)
    z = m.add_vars(2, cost=[-1, -1], ub=1)
    m.set_bounds(z, 0, [1, 0])  # unsynced columns
    m.set_cost(z[1], -5)
    assert m.solve().objective == pytest.approx(5)
    with pytest.raises(IndexError):
        m.set_row_bounds(99, 0, 1)
    with pytest.raises(IndexError):
        m.set_bounds(99, 0, 1)


def test_set_objective_after_solve_replaces_all_costs():
    m = rs.Model()
    x = m.add_vars(2, lb=1, ub=2, cost=[5, 5])
    assert m.solve().objective == pytest.approx(10)
    y = m.add_var(lb=1, ub=2, cost=7)
    m.minimize(x[0] + 100)
    assert m.solve().objective == pytest.approx(101)
    assert m.num_vars == 3 and y.index == 2


def test_solve_options():
    m, x = lp3()
    s = m.solve(time_limit=10, mip_gap=1e-6, threads=1, presolve=False, seed=3,
                options={"simplex_strategy": 1, "primal_feasibility_tolerance": 1e-8})
    assert s.status == "Optimal"
    assert m.solve().objective == pytest.approx(s.objective)
    assert m.solve(presolve="choose", options={"solver": "simplex"}).status == "Optimal"
    assert m.solve(presolve=True).status == "Optimal"
    with pytest.raises(ValueError, match="presolve"):
        m.solve(presolve="sometimes")
    with pytest.raises(TypeError, match="presolve"):
        m.solve(presolve=1.5)
    # node limit through a raw option on a MIP that needs branching
    h = hard_mip()
    assert h.solve(options={"mip_max_nodes": 1}).status in ("NodeLimit", "Optimal", "Infeasible")


def test_sync_without_solve():
    m, x = lp3()
    m._sync()
    assert m.num_constrs == 1 and m.num_nonzeros == 3
    m.add_constr(x[0] >= 1)
    m._sync()
    m._sync()
    assert m.num_constrs == 2 and m.num_nonzeros == 4
    assert m.solve().objective == pytest.approx(2)


def hard_mip(n=40, seed=1):
    """A market-split-like equality knapsack system: hard for branch and bound."""
    rng = np.random.default_rng(seed)
    m = rs.Model()
    x = m.add_vars(n, binary=True)
    for _ in range(4):
        a = rng.integers(1, 100, n)
        m.add_constr(x.dot(a.astype(float)) == float(a.sum() // 2))
    return m


def test_gil_released_during_solve():
    m = hard_mip()
    ticks = 0
    stop = threading.Event()

    def spin():
        nonlocal ticks
        while not stop.is_set():
            ticks += 1
            time.sleep(0.001)

    t = threading.Thread(target=spin)
    t.start()
    try:
        before = ticks
        sol = m.solve(time_limit=0.5, threads=1)
        during = ticks - before
    finally:
        stop.set()
        t.join()
    assert sol.time >= 0.3, sol
    assert sol.status in ("TimeLimit", "Optimal", "Infeasible")
    assert during >= 50, during  # the thread kept running during the solve


# ---------------------------------------------------------------- names


def test_names():
    m = rs.Model(name="named")
    x = m.add_vars(3, names="x")
    y = m.add_vars(2, names=["a", "b"])
    z = m.add_var(name="z")
    u = m.add_var()
    assert [m.var_name(v) for v in x] == ["x[0]", "x[1]", "x[2]"]
    assert m.var_name(y[1]) == "b" and m.var_name(z) == "z"
    r0 = m.add_constr(x.sum() <= 3, name="cap")
    rr = m.add_constrs([x[0] <= 1, x[1] <= 1], names="lim")
    r3 = m.add_constr(x[2] <= 1)
    assert m.constr_name(r0) == "cap"
    assert [m.constr_name(r) for r in rr] == ["lim[0]", "lim[1]"]
    m.solve()
    # after the solve the model has names: the unnamed get HiGHS's defaults
    assert m.var_name(x[0]) == "x[0]" and m.var_name(u) == f"c{u.index}"
    assert m.constr_name(r0) == "cap" and m.constr_name(r3) == f"r{r3}"
    w = m.add_var(name="w")
    r4 = m.add_constr(w <= 1, name="wcap")
    m.solve()
    assert m.var_name(w) == "w" and m.constr_name(r4) == "wcap"
    with pytest.raises(IndexError):
        m.constr_name(99)


def test_no_names_means_none():
    m = rs.Model()
    x = m.add_vars(2)
    m.add_constr(x[0] <= 1)
    assert m.var_name(x[0]) is None and m.constr_name(0) is None
    m.solve()
    assert m.var_name(x[0]) is None and m.constr_name(0) is None


def test_model_repr_and_counts():
    m, x = lp3()
    assert (m.num_vars, m.num_constrs, m.num_nonzeros) == (3, 1, 3)
    m.solve()
    m.add_constr(x[0] + x[1] <= 2)
    assert (m.num_vars, m.num_constrs, m.num_nonzeros) == (3, 2, 5)
    assert repr(m).startswith("Model(3 vars, 2 constrs")
