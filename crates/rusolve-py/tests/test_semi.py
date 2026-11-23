"""Semi-continuous and semi-integer columns through the Python API (`Model.set_semi`, and MPS `SC` bounds read by `Highs.readModel`). Each model is built so that
dropping the semi constraint, or treating it as a plain bound, gives a different optimum."""

import pytest

import rusolve as rs


def semi_floor():
    """min x, x semi-continuous in [3, 5]: 0 is allowed. As a plain bound the optimum is 3."""
    m = rs.Model()
    x = m.add_var(lb=3, ub=5, cost=1.0)
    return m, x


def semi_cover(integer):
    """min x + 10 y,  x + y >= 3.5,  0 <= y <= 1, x semi in [2, 6].
    Semi-continuous: x = 3.5, objective 3.5. Semi-integer: x = 4, objective 4.
    With the semi constraint dropped (x in [0, 6]) the continuous optimum is also 3.5, so the
    integer case is the one that tells semi-integer apart."""
    m = rs.Model()
    x = m.add_var(lb=2, ub=6, cost=1.0)
    y = m.add_var(lb=0, ub=1, cost=10.0)
    m.add_constr(x + y >= 3.5)
    return m, x, y


def semi_gap():
    """min x + 10 y,  x + y >= 1,  0 <= y <= 1, x semi-continuous in [3, 5].
    x = 3 costs 3, x = 0 forces y = 1 (cost 10). Dropping the lower bound gives x = 1, cost 1."""
    m = rs.Model()
    x = m.add_var(lb=3, ub=5, cost=1.0)
    y = m.add_var(lb=0, ub=1, cost=10.0)
    m.add_constr(x + y >= 1)
    return m, x, y


def test_semi_continuous_may_be_zero():
    m, x = semi_floor()
    m.set_semi(x)
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(0.0, abs=1e-9)
    assert sol.x[0] == pytest.approx(0.0, abs=1e-9)


def test_semi_continuous_keeps_its_lower_bound_when_on():
    m, x, y = semi_gap()
    m.set_semi(x)
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(3.0, abs=1e-7)
    assert sol.x[0] == pytest.approx(3.0, abs=1e-7)


def test_semi_integer_takes_an_integer_value():
    m, x, y = semi_cover(integer=True)
    m.set_semi(x, integer=True)
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(4.0, abs=1e-7)
    assert sol.x[0] == pytest.approx(4.0, abs=1e-7)


def test_semi_continuous_cover():
    m, x, y = semi_cover(integer=False)
    m.set_semi(x)
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(3.5, abs=1e-7)


def test_set_semi_after_a_solve_applies_to_the_next_solve():
    """The synced path: the column already lives in the solver when `set_semi` is called."""
    m, x = semi_floor()
    first = m.solve()
    assert first.objective == pytest.approx(3.0, abs=1e-9)
    m.set_semi(x)
    second = m.solve()
    assert second.status == "Optimal"
    assert second.objective == pytest.approx(0.0, abs=1e-9)


def test_semi_on_several_vars_at_once():
    m = rs.Model()
    xs = m.add_vars(3, lb=2, ub=4, cost=1.0)
    m.add_constr(xs[0] + xs[1] + xs[2] >= 2.5)
    m.set_semi(xs)
    sol = m.solve()
    # one column on at 2.5 is cheapest; the others stay at 0
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(2.5, abs=1e-7)
    assert sorted(round(v, 7) for v in sol.x) == [0.0, 0.0, 2.5]


MPS_SC = """NAME          SEMIGAP
ROWS
 N  obj
 G  c1
COLUMNS
    x  obj  1  c1  1
    y  obj  10  c1  1
RHS
    rhs  c1  1
BOUNDS
 LO bnd x 3
 SC bnd x 5
 UP bnd y 1
ENDATA
"""


def test_read_mps_semi_continuous(tmp_path):
    """MPS `SC` bound through the highspy-compatible reader: same model as `semi_gap`."""
    p = tmp_path / "semigap.mps"
    p.write_text(MPS_SC)
    h = rs.Highs()
    h.setOptionValue("output_flag", False)
    assert h.readModel(str(p)) == rs.HighsStatus.kOk
    assert h.run() == rs.HighsStatus.kOk
    assert h.getModelStatus() == rs.HighsModelStatus.kOptimal
    assert h.getInfo().objective_function_value == pytest.approx(3.0, abs=1e-7)
    assert h.getSolution().col_value[0] == pytest.approx(3.0, abs=1e-7)
