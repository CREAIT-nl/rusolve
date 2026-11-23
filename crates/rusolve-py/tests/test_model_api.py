"""Acceptance tests for `rusolve.Model`, walking docs/design/python-api.md section by section.

Written spec-first: where the spec leaves a detail open the test accepts either reading and
says so in a comment. Every model is small with a known optimum.
"""
import itertools

import numpy as np
import pytest
import rusolve as rs
import scipy.sparse as sp

INF = rs.inf


def solve(m, **kw):
    return m.solve(**kw)


def optimal(m, **kw):
    sol = m.solve(**kw)
    assert sol.status == "Optimal", sol.status
    assert sol.ok
    return sol


# ------------------------------------------------------------------------------ add_vars
class TestAddVars:
    def test_defaults(self):
        m = rs.Model()
        x = m.add_vars(3)
        assert m.num_vars == 3 and m.num_constrs == 0
        m.minimize(x.sum())
        sol = optimal(m)
        assert sol.objective == pytest.approx(0.0)
        assert list(sol.x) == [0.0, 0.0, 0.0]  # lb=0 default

    def test_default_ub_is_inf(self):
        m = rs.Model()
        x = m.add_vars(1)
        m.maximize(x[0])
        assert solve(m).status == "Unbounded"

    def test_scalar_lb_ub_cost(self):
        m = rs.Model()
        m.add_vars(2, lb=1, ub=5, cost=2)
        sol = optimal(m)  # minimise 2 x0 + 2 x1 -> both at lb
        assert sol.objective == pytest.approx(4.0)

    def test_sequence_bounds_and_cost(self):
        m = rs.Model()
        m.add_vars(2, lb=[1, 2], ub=[3, 4], cost=[1, -1])
        sol = optimal(m)  # x0 -> 1, x1 -> 4
        assert list(sol.x) == pytest.approx([1.0, 4.0])
        assert sol.objective == pytest.approx(-3.0)

    @pytest.mark.parametrize("dtype", [np.float64, np.float32, np.int64])
    def test_numpy_bounds_and_cost(self, dtype):
        m = rs.Model()
        m.add_vars(2, lb=np.array([1, 2], dtype=dtype), ub=np.array([3, 4], dtype=dtype),
                   cost=np.array([1, -1], dtype=dtype))
        sol = optimal(m)
        assert list(sol.x) == pytest.approx([1.0, 4.0])

    def test_infinite_lower_bound(self):
        m = rs.Model()
        x = m.add_vars(1, lb=-INF, ub=2, cost=1)
        m.add_constr(x[0] >= -7)
        sol = optimal(m)
        assert sol.x[0] == pytest.approx(-7.0)

    def test_integer_flag(self):
        m = rs.Model()
        x = m.add_vars(1, ub=2.5, integer=True)
        m.maximize(x[0])
        sol = optimal(m)
        assert sol.x[0] == pytest.approx(2.0)

    def test_integer_sequence(self):
        m = rs.Model()
        x = m.add_vars(2, ub=2.5, integer=[True, False])
        m.maximize(x[0] + x[1])
        sol = optimal(m)
        assert list(sol.x) == pytest.approx([2.0, 2.5])

    def test_binary(self):
        m = rs.Model()
        x = m.add_vars(3, binary=True)
        m.add_constr(x.sum() <= 2)
        m.maximize(x.sum())
        sol = optimal(m)
        assert sol.objective == pytest.approx(2.0)
        assert all(v in (0.0, 1.0) for v in sol.x)

    def test_binary_is_integer_with_default_bounds_01(self):
        m = rs.Model()
        x = m.add_vars(1, binary=True, ub=0.5, cost=-1)  # an explicit bound still applies
        assert optimal(m)[x[0]] == pytest.approx(0.0)

    def test_add_var(self):
        m = rs.Model()
        a = m.add_var(lb=1, ub=4, cost=3)
        b = m.add_var(ub=2, cost=-1, integer=True)
        assert (int(a), int(b)) == (0, 1)
        assert hash(a) == hash(0) and hash(b) == hash(1)
        sol = optimal(m)
        assert list(sol.x) == pytest.approx([1.0, 2.0])
        assert sol.objective == pytest.approx(1.0)

    def test_add_var_binary(self):
        m = rs.Model()
        v = m.add_var(binary=True, cost=-1)
        assert optimal(m)[v] == pytest.approx(1.0)

    def test_names(self):
        m = rs.Model()
        a = m.add_vars(2, names=["u", "v"])
        b = m.add_vars(2, names="p")
        c = m.add_var(name="solo")
        assert m.var_name(a[1]) == "v"
        assert m.var_name(b[0]) == "p[0]" and m.var_name(b[1]) == "p[1]"
        assert m.var_name(c) == "solo"

    def test_default_names_follow_highs_prefixes(self):
        m = rs.Model()
        m.add_vars(1)
        y = m.add_var(name="y")
        z = m.add_var()
        m.add_constr({z: 1}, "<=", 1)
        r = m.add_constr({z: 1}, ">=", 0, name="named")
        assert m.var_name(y) == "y" and m.constr_name(r) == "named"
        m.solve()  # "unnamed columns and rows get c<j> / r<i> at build time"
        assert m.var_name(z) == "c2"
        assert m.constr_name(0) == "r0"

    def test_length_mismatch(self):
        m = rs.Model()
        with pytest.raises((ValueError, TypeError)):
            m.add_vars(3, lb=[0, 0])


class TestVarsBlock:
    def test_range_and_indexing(self):
        m = rs.Model()
        m.add_vars(2)
        x = m.add_vars(4)
        assert x.start == 2 and x.stop == 6 and len(x) == 4
        assert list(x.indices) == [2, 3, 4, 5]
        assert [int(v) for v in x] == [2, 3, 4, 5]
        assert int(x[1]) == 3 and int(x[-1]) == 5
        sub = x[1:3]
        assert isinstance(sub, type(x)) and sub.start == 3 and sub.stop == 5

    def test_sum_and_dot(self):
        m = rs.Model()
        x = m.add_vars(3, ub=1)
        m.maximize(x.dot([1, 2, 3]))
        m.add_constr(x.sum() <= 2)
        sol = optimal(m)
        assert sol.objective == pytest.approx(5.0)
        m.maximize(x.dot(np.array([3.0, 2.0, 1.0])))
        assert optimal(m).objective == pytest.approx(5.0)


# ------------------------------------------------------------------------- constraints
def two_var_model():
    """max x + y over [0, inf)^2, solved under a row the test adds."""
    m = rs.Model()
    x = m.add_vars(2)
    m.maximize(x[0] + x[1])
    return m, x


class TestRows:
    def test_dict_row_var_keys(self):
        m, x = two_var_model()
        r = m.add_constr({x[0]: 1.0, x[1]: 2.0}, "<=", 4.0)
        assert r == 0 and isinstance(r, int)
        m.add_constr({x[0]: 1.0}, "<=", 3.0)
        assert optimal(m).objective == pytest.approx(3.5)

    def test_dict_row_int_keys(self):
        m, x = two_var_model()
        m.add_constr({0: 1.0, 1: 2.0}, "<=", 4.0)
        m.add_constr({0: 1}, "<=", 3)
        assert optimal(m).objective == pytest.approx(3.5)

    def test_index_value_pair(self):
        m, x = two_var_model()
        m.add_constr(([0, 1], [1.0, 2.0]), "<=", 4.0)
        m.add_constr((np.array([0], dtype=np.int32), np.array([1.0])), "<=", 3.0)
        m.add_constr((np.array([1], dtype=np.int64), np.array([1.0])), "<=", 10.0)
        assert optimal(m).objective == pytest.approx(3.5)

    def test_linexpr_and_var_rows(self):
        m, x = two_var_model()
        m.add_constr(x[0] + 2 * x[1], "<=", 4.0)
        m.add_constr(x[0], "<=", 3.0)
        assert optimal(m).objective == pytest.approx(3.5)

    @pytest.mark.parametrize("sense,rhs,expect", [
        ("<=", 4, 4), ("L", 4, 4), ("==", 4, 4), ("E", 4, 4)])
    def test_sense_spellings_upper(self, sense, rhs, expect):
        m, x = two_var_model()
        m.add_constr({x[0]: 1, x[1]: 1}, sense, rhs)
        assert optimal(m).objective == pytest.approx(expect)

    @pytest.mark.parametrize("sense", [">=", "G"])
    def test_sense_ge(self, sense):
        m = rs.Model()
        x = m.add_vars(2, cost=[1, 2])
        m.add_constr({x[0]: 1, x[1]: 1}, sense, 3)
        assert optimal(m).objective == pytest.approx(3.0)

    def test_equality(self):
        m = rs.Model()
        x = m.add_vars(2, cost=[1, 2])
        m.add_constr(x[0] + x[1] == 3)
        m.add_constr(x[0] <= 1)
        sol = optimal(m)
        assert list(sol.x) == pytest.approx([1.0, 2.0])

    def test_constraint_object(self):
        m, x = two_var_model()
        r = m.add_constr(x[0] + 2 * x[1] <= 4)
        m.add_constr(x[0] <= 3)
        assert r == 0
        assert optimal(m).objective == pytest.approx(3.5)

    def test_constraint_constant_moves_to_bounds(self):
        m, x = two_var_model()
        m.add_constr(x[0] + 1 <= 3)  # x0 <= 2
        m.add_constr(x[1] <= 1)
        assert optimal(m).objective == pytest.approx(3.0)

    def test_reversed_comparison(self):
        m, x = two_var_model()
        m.add_constr(5 >= x[0] + x[1])
        assert optimal(m).objective == pytest.approx(5.0)

    def test_ranged_lo_hi(self):
        m = rs.Model()
        x = m.add_vars(1, cost=1)
        m.add_constr({x[0]: 1}, lo=1.0, hi=4.0)
        assert optimal(m).objective == pytest.approx(1.0)
        m.maximize(x[0])
        assert optimal(m).objective == pytest.approx(4.0)

    def test_between(self):
        m = rs.Model()
        x = m.add_vars(2, cost=[1, 1])
        m.add_constr((x[0] + x[1]).between(2, 6))
        assert optimal(m).objective == pytest.approx(2.0)
        m.maximize(x[0] + x[1])
        assert optimal(m).objective == pytest.approx(6.0)

    def test_duplicates_summed_and_zeros_dropped(self):
        m, x = two_var_model()
        # x0 + x0 + x1 + 0*x1: a 2 x0 + x1 row
        m.add_constr(x[0] + x[0] + x[1] + 0 * x[1], "<=", 4)
        m.add_constr(x[0] <= 1)
        m.maximize(3 * x[0] + x[1])
        assert optimal(m).objective == pytest.approx(5.0)  # x0=1, x1=2 (unmerged: 6)
        m.add_constr(([0, 0], [1.0, 1.0]), "<=", 1.0)  # 2 x0 <= 1
        assert optimal(m).objective == pytest.approx(1.5 + 3.0)

    def test_row_names(self):
        m, x = two_var_model()
        r = m.add_constr(x[0] + x[1] <= 1, name="cap")
        assert m.constr_name(r) == "cap"
        assert m.num_constrs == 1 and m.num_nonzeros == 2


class TestAddConstrs:
    def test_scalar_rhs_and_sense(self):
        m = rs.Model()
        x = m.add_vars(3)
        m.maximize(x.sum())
        rr = m.add_constrs([{x[j]: 1} for j in range(3)], "<=", 2.0)
        assert len(rr) == 3 and list(rr) == [0, 1, 2]
        assert optimal(m).objective == pytest.approx(6.0)

    def test_vector_rhs(self):
        m = rs.Model()
        x = m.add_vars(3)
        m.maximize(x.sum())
        m.add_constrs([{x[j]: 1} for j in range(3)], "<=", [1, 2, 3])
        assert optimal(m).objective == pytest.approx(6.0)

    def test_numpy_rhs(self):
        m = rs.Model()
        x = m.add_vars(3)
        m.maximize(x.sum())
        m.add_constrs([{x[j]: 1} for j in range(3)], "<=", np.array([1.0, 2.0, 3.0]))
        assert optimal(m).objective == pytest.approx(6.0)

    def test_vector_sense(self):
        m = rs.Model()
        x = m.add_vars(3, cost=1)
        rows = [{x[j]: 1} for j in range(3)]
        m.add_constrs(rows, ["==", ">=", "<="], [1, 2, 3])
        sol = optimal(m)
        assert list(sol.x) == pytest.approx([1.0, 2.0, 0.0])

    def test_constraint_objects(self):
        m = rs.Model()
        x = m.add_vars(3)
        m.maximize(x.sum())
        rr = m.add_constrs(x[j] <= j + 1 for j in range(3))
        assert len(rr) == 3
        assert optimal(m).objective == pytest.approx(6.0)

    def test_names(self):
        m = rs.Model()
        x = m.add_vars(2)
        rr = m.add_constrs([{x[0]: 1}, {x[1]: 1}], "<=", 1, names=["a", "b"])
        assert [m.constr_name(r) for r in rr] == ["a", "b"]

    def test_range_continues_numbering(self):
        m = rs.Model()
        x = m.add_vars(2)
        m.add_constr({x[0]: 1}, "<=", 1)
        rr = m.add_constrs([{x[1]: 1}, {x[0]: 1, x[1]: 1}], "<=", 1)
        assert (rr.start, rr.stop) == (1, 3)

    def test_dict_rows_with_int_keys_and_empty_iterable(self):
        m = rs.Model()
        m.add_vars(2)
        rr = m.add_constrs([], "<=", [])
        assert len(rr) == 0
        rr = m.add_constrs([{0: 1.0, 1: 1.0}], "<=", 1.0)
        assert len(rr) == 1


def small_matrix():
    # rows: x0 + 2 x1 <= 4 ; 3 x0 + x1 <= 6 ; max x0 + x1 -> (1.6, 1.2), obj 2.8
    A = np.array([[1.0, 2.0], [3.0, 1.0]])
    return A, np.array([4.0, 6.0])


class TestAddMatrix:
    @pytest.mark.parametrize("fmt", ["csr", "csc", "coo", "dense", "csr_array", "list"])
    def test_formats(self, fmt):
        A, b = small_matrix()
        M = {"csr": sp.csr_matrix(A), "csc": sp.csc_matrix(A), "coo": sp.coo_matrix(A),
             "dense": A, "csr_array": sp.csr_array(A), "list": None}[fmt]
        if fmt == "list":
            pytest.skip("lists are not a spec'd matrix form (2-D float64 buffer / scipy.sparse)")
        m = rs.Model()
        x = m.add_vars(2)
        m.maximize(x.sum())
        rr = m.add_matrix(M, "<=", b)
        assert len(rr) == 2
        assert m.num_nonzeros == 4
        sol = optimal(m)
        assert sol.objective == pytest.approx(2.8)
        assert list(sol.x) == pytest.approx([1.6, 1.2])

    def test_dense_keeps_nonzeros_only(self):
        A = np.array([[1.0, 0.0, 2.0], [0.0, 0.0, 1.0]])
        m = rs.Model()
        m.add_vars(3)
        m.add_matrix(A, "<=", [1.0, 1.0])
        assert m.num_nonzeros == 3

    def test_sparse_with_explicit_zeros(self):
        A = sp.csr_matrix(
            (np.array([1.0, 0.0, 2.0]), np.array([0, 1, 1]), np.array([0, 2, 3])), shape=(2, 2))
        m = rs.Model()
        x = m.add_vars(2)
        m.maximize(x.sum())
        m.add_matrix(A, "<=", [1.0, 2.0])
        sol = optimal(m)
        # explicit zero never constrains: x0 <= 1, 2 x1 <= 2
        assert sol.objective == pytest.approx(2.0)

    def test_senses_and_rhs_vector(self):
        A, b = small_matrix()
        m = rs.Model()
        m.add_vars(2, cost=1)
        m.add_matrix(sp.csr_matrix(A), ">=", [1.0, 3.0])
        sol = optimal(m)
        # x0 + 2 x1 >= 1, 3 x0 + x1 >= 3 -> min x0 + x1: vertex (1, 0) gives 1
        assert sol.objective == pytest.approx(1.0)

    def test_vars_mapping(self):
        A = np.array([[1.0, 1.0]])
        m = rs.Model()
        m.add_vars(1)  # column 0 is not part of the matrix
        y = m.add_vars(2)
        m.maximize(y.sum())
        m.add_matrix(A, "<=", [3.0], vars=y)
        sol = optimal(m)
        assert sol.objective == pytest.approx(3.0)
        m2 = rs.Model()
        v = m2.add_vars(3)
        m2.maximize(v.sum())
        m2.add_matrix(A, "<=", [3.0], vars=[v[2], v[0]])
        m2.add_constr(v[1] <= 1)
        assert optimal(m2).objective == pytest.approx(4.0)

    def test_names(self):
        A, b = small_matrix()
        m = rs.Model()
        m.add_vars(2)
        rr = m.add_matrix(A, "<=", b, names=["a", "b"])
        assert [m.constr_name(r) for r in rr] == ["a", "b"]

    def test_wrong_shape(self):
        m = rs.Model()
        m.add_vars(2)
        with pytest.raises((ValueError, IndexError)):
            m.add_matrix(np.ones((2, 5)), "<=", [1.0, 1.0])


# -------------------------------------------------------------------------- objective
class TestObjective:
    def setup_method(self):
        self.m = rs.Model()
        self.x = self.m.add_vars(3, ub=2)

    def test_dict(self):
        x = self.x
        self.m.minimize({x[0]: -1, x[2]: -2})
        assert optimal(self.m).objective == pytest.approx(-6.0)

    def test_dict_int_keys(self):
        self.m.maximize({0: 1, 1: 1})
        assert optimal(self.m).objective == pytest.approx(4.0)

    def test_linexpr(self):
        x = self.x
        self.m.maximize(x[0] + 2 * x[1] - x[2])
        assert optimal(self.m).objective == pytest.approx(6.0)

    def test_var(self):
        self.m.maximize(self.x[1])
        assert optimal(self.m).objective == pytest.approx(2.0)

    def test_pair(self):
        self.m.maximize(([0, 2], [1.0, 3.0]))
        assert optimal(self.m).objective == pytest.approx(8.0)

    def test_dense_sequence_and_buffer(self):
        self.m.maximize([1.0, 2.0, 3.0])
        assert optimal(self.m).objective == pytest.approx(12.0)
        self.m.minimize(np.array([-1.0, 0.0, -1.0]))
        assert optimal(self.m).objective == pytest.approx(-4.0)

    def test_constant_becomes_offset(self):
        x = self.x
        self.m.maximize(x[0] + 10)
        assert optimal(self.m).objective == pytest.approx(12.0)

    def test_replaces_whole_objective(self):
        m = rs.Model()
        y = m.add_vars(2, ub=1, cost=[5, 5])
        m.maximize(y[0])  # y1's cost of 5 is cleared
        assert optimal(m).objective == pytest.approx(1.0)

    def test_set_objective_sense(self):
        x = self.x
        self.m.set_objective(x[0], sense="max")
        assert optimal(self.m).objective == pytest.approx(2.0)
        self.m.set_objective(x[0], sense="min")
        assert optimal(self.m).objective == pytest.approx(0.0)

    def test_model_sense_argument(self):
        m = rs.Model(sense="max", name="demo")
        m.add_vars(1, ub=3, cost=1)
        assert optimal(m).objective == pytest.approx(3.0)

    def test_cost_in_add_vars_is_the_objective(self):
        m = rs.Model()
        m.add_vars(2, ub=1, cost=[1, 2])
        assert optimal(m).objective == pytest.approx(0.0)


# ----------------------------------------------------------------------------- solution
def lp_with_duals():
    # min 2x + 3y  s.t. x + y >= 4, x <= 3  ->  x=3, y=1, obj 9; dual of row 3, rc(x) = -1
    m = rs.Model()
    x = m.add_vars(2, cost=[2, 3], ub=[3, INF])
    m.add_constr(x[0] + x[1] >= 4)
    return m, x


class TestSolution:
    def test_fields_lp(self):
        m, x = lp_with_duals()
        sol = optimal(m)
        assert sol.status == "Optimal" and sol.ok is True
        assert sol.objective == pytest.approx(9.0)
        assert list(sol.x) == pytest.approx([3.0, 1.0])
        assert list(sol.row_activity) == pytest.approx([4.0])
        assert list(sol.duals) == pytest.approx([3.0])
        assert list(sol.reduced_costs) == pytest.approx([-1.0, 0.0])
        assert sol.iterations >= 0 and sol.time >= 0.0
        assert sol.nodes in (0, None) or sol.nodes >= 0

    def test_x_is_numpy_float64(self):
        m, _ = lp_with_duals()
        x = optimal(m).x
        assert isinstance(x, np.ndarray) and x.dtype == np.float64 and x.shape == (2,)

    def test_mip_fields(self):
        m = rs.Model()
        x = m.add_vars(3, binary=True, cost=[-4, -5, -3])
        m.add_constr(2 * x[0] + 3 * x[1] + 2 * x[2] <= 4)
        sol = optimal(m)
        assert sol.objective == pytest.approx(-7.0)  # x0 + x2 -> 7 ; x1 alone -> 5
        assert sol.gap is None or sol.gap <= 1e-4
        assert sol.bound is None or sol.bound <= sol.objective + 1e-6
        assert sol.nodes is None or sol.nodes >= 0
        assert sol.duals is None  # duals: LP only

    def test_indexing(self):
        m, x = lp_with_duals()
        sol = optimal(m)
        assert sol[x[0]] == pytest.approx(3.0)
        assert sol[1] == pytest.approx(1.0)
        assert sol[2 * x[0] + x[1]] == pytest.approx(7.0)
        assert sol[{x[0]: 1.0, x[1]: 10.0}] == pytest.approx(13.0)

    def test_expr_constant_is_included(self):
        # the spec lists LinExpr as an index; a constant is part of the expression's value
        m, x = lp_with_duals()
        sol = optimal(m)
        assert sol[x[0] + 5] == pytest.approx(8.0)

    def test_values_helper(self):
        m = rs.Model()
        x = m.add_vars(3, ub=[1, 2, 3], cost=-1)
        sol = optimal(m)
        mapping = {"a": x[0], "b": x[1], "c": x[2]}
        assert sol.values(mapping) == {"a": pytest.approx(1.0), "b": pytest.approx(2.0),
                                       "c": pytest.approx(3.0)}
        assert set(sol.values(mapping, above=1.5)) == {"b", "c"}

    def test_selected_helper(self):
        m = rs.Model()
        x = m.add_vars(4, binary=True, cost=[-1, 1, -2, 1])
        sol = optimal(m)
        mapping = {f"i{j}": x[j] for j in range(4)}
        assert sorted(sol.selected(mapping)) == ["i0", "i2"]
        assert sorted(sol.selected(mapping, threshold=-1)) == ["i0", "i1", "i2", "i3"]

    def test_infeasible(self):
        m = rs.Model()
        x = m.add_vars(1)
        m.add_constr(x[0] >= 2)
        m.add_constr(x[0] <= 1)
        sol = solve(m)
        assert sol.status == "Infeasible"
        assert sol.ok is False and sol.x is None

    def test_unbounded(self):
        m = rs.Model()
        x = m.add_vars(2)
        m.add_constr(x[0] - x[1] <= 1)
        m.maximize(x[0] + x[1])
        sol = solve(m)
        assert sol.status == "Unbounded"
        assert sol.ok is False

    def test_infeasible_mip(self):
        m = rs.Model()
        x = m.add_vars(1, ub=1, integer=True)
        m.add_constr(2 * x[0] == 1)
        assert solve(m).status == "Infeasible"

    def test_time_limit_option_accepted(self):
        m, _ = lp_with_duals()
        sol = m.solve(time_limit=10.0, mip_gap=0.0, threads=1, presolve=True, seed=1)
        assert sol.status == "Optimal"

    def test_threads_can_change_between_solves(self):
        # The executor is per OS thread; solving with 1, 4, 2 and 1
        # threads in one process used to fail with "global scheduler has already been
        # initialized to use 1 threads".
        m = rs.Model()
        x = m.add_vars(3, ub=4, integer=True)
        m.add_constr(x[0] + 2 * x[1] + x[2] <= 8)
        m.add_constr(x[2] + 3 * x[1] <= 10)
        m.maximize(3 * x[0] + 4 * x[1] - 5 * x[2])
        objs = []
        for threads in (1, 4, 2, 1, 4, 1):
            sol = m.solve(threads=threads)
            assert sol.status == "Optimal", (threads, sol.status)
            objs.append(sol.objective)
        assert max(objs) - min(objs) < 1e-9, objs

    def test_raw_options(self):
        m, _ = lp_with_duals()
        sol = m.solve(options={"presolve": "off"})
        assert sol.status == "Optimal"
        with pytest.raises((ValueError, KeyError, TypeError, RuntimeError)):
            m.solve(options={"no_such_option": 1})

    def test_verbose_prints_log(self, capfd):
        m, _ = lp_with_duals()
        m.solve(verbose=True)
        out = capfd.readouterr()
        assert out.out or out.err

    def test_default_is_quiet(self, capfd):
        m, _ = lp_with_duals()
        m.solve()
        out = capfd.readouterr()
        assert not out.out and not out.err

    def test_mip_against_brute_force(self):
        values = [10, 13, 18, 31, 7, 15]
        weights = [11, 15, 20, 35, 10, 33]
        cap = 47
        best = max(sum(v for v, s in zip(values, pick) if s)
                   for pick in itertools.product([0, 1], repeat=6)
                   if sum(w for w, s in zip(weights, pick) if s) <= cap)
        m = rs.Model(sense="max")
        x = m.add_vars(6, binary=True, cost=values)
        m.add_constr(x.dot(weights) <= cap)
        assert optimal(m).objective == pytest.approx(best)


# ----------------------------------------------------------------------------- changes
class TestChanges:
    def test_bounds_and_resolve(self):
        m = rs.Model()
        x = m.add_vars(2, ub=1, cost=[-1, -1])
        m.add_constr(x[0] + x[1] <= 5)
        assert optimal(m).objective == pytest.approx(-2.0)
        m.set_bounds(x, 0, 3)
        assert optimal(m).objective == pytest.approx(-5.0)
        m.set_bounds([x[0]], [0.0], [0.5])
        assert optimal(m).objective == pytest.approx(-3.5)

    def test_set_cost(self):
        m = rs.Model()
        x = m.add_vars(2, ub=1, cost=[-1, -1])
        assert optimal(m).objective == pytest.approx(-2.0)
        m.set_cost(x, [1, -3])
        sol = optimal(m)
        assert sol.objective == pytest.approx(-3.0)
        assert list(sol.x) == pytest.approx([0.0, 1.0])

    def test_set_integer(self):
        m = rs.Model()
        x = m.add_vars(1, ub=2.5)
        m.maximize(x[0])
        assert optimal(m).objective == pytest.approx(2.5)
        m.set_integer(x, True)
        assert optimal(m).objective == pytest.approx(2.0)
        m.set_integer(x, False)
        assert optimal(m).objective == pytest.approx(2.5)

    def test_set_row_bounds(self):
        m = rs.Model()
        x = m.add_vars(1)
        m.maximize(x[0])
        r = m.add_constr(x[0] <= 4)
        assert optimal(m).objective == pytest.approx(4.0)
        m.set_row_bounds([r], [-INF], [7.0])
        assert optimal(m).objective == pytest.approx(7.0)
        m.set_row_bounds([r], [1.0], [2.0])
        assert optimal(m).objective == pytest.approx(2.0)

    def test_set_coef(self):
        m = rs.Model()
        x = m.add_vars(2)
        m.maximize(x[0] + x[1])
        r = m.add_constr(x[0] + x[1] <= 4)
        m.add_constr(x[0] <= 1)
        assert optimal(m).objective == pytest.approx(4.0)
        m.set_coef(r, x[1], 2.0)  # x0 + 2 x1 <= 4
        assert optimal(m).objective == pytest.approx(2.5)

    def test_add_after_solve(self):
        m = rs.Model()
        x = m.add_vars(2, ub=10)
        m.maximize(x[0] + x[1])
        m.add_constr(x[0] + x[1] <= 8)
        assert optimal(m).objective == pytest.approx(8.0)
        y = m.add_var(ub=1, cost=1)  # new column after a solve, objective via cost
        m.add_constr(x[0] <= 3)
        m.add_constr(x[0] + x[1] + y <= 8.5)
        m.maximize(x[0] + x[1] + y)
        sol = optimal(m)
        assert sol.objective == pytest.approx(8.5)
        assert m.num_vars == 3 and m.num_constrs == 3
        assert sol[y] <= 1.0 + 1e-9

    def test_change_objective_after_solve(self):
        m = rs.Model()
        x = m.add_vars(2, ub=1)
        m.maximize(x[0])
        assert optimal(m).objective == pytest.approx(1.0)
        m.minimize(x[0] + 7)
        assert optimal(m).objective == pytest.approx(7.0)

    def test_resolve_matches_fresh_model(self):
        rng = np.random.default_rng(3)
        A = rng.uniform(0, 1, (6, 5))
        b = rng.uniform(5, 9, 6)
        c = rng.uniform(1, 2, 5)

        def fresh(ub):
            m = rs.Model(sense="max")
            x = m.add_vars(5, ub=ub, cost=c)
            m.add_matrix(A, "<=", b)
            return m, x

        m, x = fresh(1.0)
        first = optimal(m).objective
        m.set_bounds(x, 0, 2.0)
        again = optimal(m).objective
        ref = optimal(fresh(2.0)[0]).objective
        assert again == pytest.approx(ref)
        assert again >= first - 1e-9

    def test_counts(self):
        m = rs.Model()
        x = m.add_vars(3)
        m.add_constr(x[0] + x[1], "<=", 1)
        m.add_constr({x[2]: 4}, ">=", 0)
        assert (m.num_vars, m.num_constrs, m.num_nonzeros) == (3, 2, 3)


# -------------------------------------------------------------------------------- errors
def build_and_solve(row_factory):
    m = rs.Model()
    x = m.add_vars(2)
    n0 = m.num_constrs
    try:
        row_factory(m, x)
    finally:
        # errors are raised at the add_* call, and a failed call adds nothing
        assert m.num_constrs == n0


class TestErrors:
    def test_nan_coefficient_names_the_row(self):
        def build(m, x):
            m.add_constr({x[0]: float("nan")}, "<=", 1.0, name="badrow")

        with pytest.raises(ValueError) as ei:
            build_and_solve(build)
        assert "badrow" in str(ei.value)

    def test_nan_rhs(self):
        def build(m, x):
            m.add_constr({x[0]: 1.0}, "<=", float("nan"))

        with pytest.raises(ValueError):
            build_and_solve(build)

    def test_nan_in_bounds_and_cost(self):
        m = rs.Model()
        with pytest.raises(ValueError):
            m.add_vars(1, lb=float("nan"))
        m = rs.Model()
        with pytest.raises(ValueError):
            m.add_vars(1, cost=float("nan"))

    def test_nan_in_objective_and_matrix(self):
        m = rs.Model()
        x = m.add_vars(2)
        with pytest.raises(ValueError):
            m.minimize({x[0]: float("nan")})
        m = rs.Model()
        m.add_vars(2)
        with pytest.raises(ValueError):
            m.add_matrix(np.array([[1.0, float("nan")]]), "<=", [1.0])

    def test_infinite_coefficient(self):
        def build(m, x):
            m.add_constr({x[0]: INF}, "<=", 1.0)

        with pytest.raises(ValueError):
            build_and_solve(build)

    def test_var_from_another_model(self):
        other = rs.Model()
        y = other.add_vars(2)
        m = rs.Model()
        m.add_vars(2)
        with pytest.raises(ValueError):
            m.add_constr(y[0] + y[1] <= 1)
        with pytest.raises(ValueError):
            m.minimize(y[0])

    def test_bad_sense(self):
        m = rs.Model()
        x = m.add_vars(1)
        with pytest.raises(ValueError):
            m.add_constr({x[0]: 1}, "<>", 1.0)

    def test_bad_model_sense(self):
        with pytest.raises(ValueError):
            rs.Model(sense="sideways")

    def test_index_out_of_range(self):
        m = rs.Model()
        m.add_vars(2)
        with pytest.raises((IndexError, ValueError)):
            m.add_constr({5: 1.0}, "<=", 1.0)

    def test_bool_of_constraint_raises(self):
        m = rs.Model()
        x = m.add_vars(2)
        with pytest.raises(Exception) as ei:
            if x[0] <= 1:
                pass
        assert "between" in str(ei.value)
        with pytest.raises(Exception):
            bool(x[0] + x[1] >= 2)

    def test_var_equality_is_index_equality(self):
        m = rs.Model()
        x = m.add_vars(2)
        assert bool(x[0] == x[0]) is True
        assert bool(x[0] == x[1]) is False
        assert x[0] in {x[0]: 1}
        assert x[1] not in {x[0]: 1}
        assert (x[0] in [x[1], x[0]]) is True

    def test_var_equals_number_is_a_constraint(self):
        m = rs.Model()
        x = m.add_vars(1)
        m.minimize(x[0])
        m.add_constr(x[0] == 3)
        assert optimal(m).x[0] == pytest.approx(3.0)

    def test_non_numeric_coefficient(self):
        m = rs.Model()
        x = m.add_vars(1)
        with pytest.raises(TypeError):
            m.add_constr({x[0]: "two"}, "<=", 1.0)

    def test_division_by_number(self):
        m = rs.Model()
        x = m.add_vars(1, ub=1)
        m.maximize(x[0] / 4)
        assert optimal(m).objective == pytest.approx(0.25)

    def test_nan_is_not_silently_ignored_in_expressions(self):
        m = rs.Model()
        x = m.add_vars(1)
        with pytest.raises(ValueError):
            m.add_constr(float("nan") * x[0] <= 1)


# ----------------------------------------------------------------------------- expression
class TestExpressions:
    def test_arithmetic(self):
        m = rs.Model()
        x = m.add_vars(3, ub=1)
        e = 2 * x[0] - x[1] + 3 * x[2] / 3 + 1 - (-x[1])
        m.maximize(e)  # 2 x0 + 0 x1 + x2 + 1
        assert optimal(m).objective == pytest.approx(4.0)

    def test_inplace_ops(self):
        m = rs.Model()
        x = m.add_vars(3, ub=1)
        e = x[0] + 0
        e += x[1]
        e += 2 * x[2]
        e -= x[0]
        e += 5
        m.maximize(e)
        assert optimal(m).objective == pytest.approx(8.0)

    def test_radd_and_sum_builtin(self):
        m = rs.Model()
        x = m.add_vars(3, ub=1)
        m.maximize(sum(x))  # builtin sum starts at int 0: needs __radd__
        assert optimal(m).objective == pytest.approx(3.0)

    def test_quicksum_vars_and_numbers(self):
        m = rs.Model()
        x = m.add_vars(3, ub=1)
        m.maximize(rs.quicksum([x[0], x[1], 2.5]))
        assert optimal(m).objective == pytest.approx(4.5)

    def test_quicksum_generator_of_products(self):
        m = rs.Model()
        x = m.add_vars(4, ub=1)
        m.maximize(rs.quicksum((j + 1) * x[j] for j in range(4)))
        assert optimal(m).objective == pytest.approx(10.0)

    def test_quicksum_coef_var_pairs(self):
        m = rs.Model()
        x = m.add_vars(3, ub=1)
        m.maximize(rs.quicksum((c, x[j]) for j, c in enumerate([1.0, 2.0, 3.0])))
        assert optimal(m).objective == pytest.approx(6.0)

    def test_quicksum_of_linexprs(self):
        m = rs.Model()
        x = m.add_vars(2, ub=1)
        m.maximize(rs.quicksum([x[0] + x[1], 2 * x[0] + 1]))
        assert optimal(m).objective == pytest.approx(5.0)

    def test_quicksum_empty(self):
        m = rs.Model()
        m.add_vars(1, ub=1)
        m.minimize(rs.quicksum([]))
        assert optimal(m).objective == pytest.approx(0.0)

    def test_constraint_from_expressions(self):
        m = rs.Model()
        x = m.add_vars(2)
        m.maximize(x[0])
        m.add_constr(x[0] <= x[1] + 1)
        m.add_constr(x[1] <= 2)
        assert optimal(m).objective == pytest.approx(3.0)

    def test_ge_and_eq_constraints(self):
        m = rs.Model()
        x = m.add_vars(2, cost=[1, 1])
        m.add_constr(x[0] + x[1] >= 2)
        m.add_constr(x[0] == 0.5)
        sol = optimal(m)
        assert list(sol.x) == pytest.approx([0.5, 1.5])


# ---------------------------------------------------------------------- end to end sample
def test_readme_example():
    m = rs.Model()
    x = m.add_vars(3, binary=True)
    m.add_constr(x[0] + x[1] + x[2] <= 2)
    m.maximize(5 * x[0] + 4 * x[1] + 3 * x[2])
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(9.0)
