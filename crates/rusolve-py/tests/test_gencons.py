"""General constraints through the Python API: indicator, SOS1/SOS2, piecewise-linear, abs/min/max
and and/or/not (`Model.add_*`). Every model is built so that dropping the object, or treating it
as a plain row, gives a different optimum. The solve goes through the same linearization as
`rusolve mip` on a file holding the same objects, and `check` and the result envelope report on
the objects, not on generated rows."""

import json
import os
import pathlib
import subprocess

import pytest

import rusolve as rs

ROOT = pathlib.Path(__file__).resolve().parents[3]
CLI = pathlib.Path(os.environ.get("CARGO_TARGET_DIR", ROOT / "target")) / "release/rusolve-cli"


def cli_objective(path):
    if not CLI.exists():
        pytest.skip("rusolve-cli is not built")
    out = subprocess.run([str(CLI), "mip", str(path), "--format", "json", "--threads", "1"],
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out)["objective"]


# ---------------------------------------------------------------- builders


def indicator_model():
    """min 2x + 7z + y, x + y >= 3, z = 0 -> x >= 4. z = 0 costs 8, z = 1 costs 10."""
    m = rs.Model(name="ind")
    z = m.add_var(binary=True, cost=7, name="z")
    x = m.add_var(ub=10, cost=2, name="x")
    y = m.add_var(ub=10, cost=1, name="y")
    m.add_constr(x + y >= 3)
    name = m.add_indicator(z, 0, x >= 4, name="need_x")
    return m, (z, x, y), name


def bigm_twin():
    m = rs.Model()
    z = m.add_var(binary=True, cost=7, name="z")
    x = m.add_var(ub=10, cost=2, name="x")
    y = m.add_var(ub=10, cost=1, name="y")
    m.add_constr(x + y >= 3)
    m.add_constr(x + 4 * z >= 4)
    return m


def sos1_model():
    """max x1 + 2 x2 + 3 x3 over [0, 1], at most one nonzero: 3 (without the set: 6)."""
    m = rs.Model(sense="max")
    x = m.add_vars(3, ub=1, names="x")
    m.maximize(x[0] + 2 * x[1] + 3 * x[2])
    m.add_sos1(x)
    return m


def sos2_model():
    """max 5 x0 + x1 + 5 x2 over [0, 1], two adjacent nonzero: 6 (without the set: 11)."""
    m = rs.Model(sense="max")
    x = m.add_vars(3, ub=1, names="x")
    m.maximize(5 * x[0] + x[1] + 5 * x[2])
    m.add_sos2(x, weights=[1, 2, 3])
    return m


def pwl_model():
    """min y, y = f(x) non-convex through (0,0) (2,4) (4,4.5) (6,10), x >= 3: f(3) = 4.25."""
    m = rs.Model()
    x = m.add_var(ub=6, name="x")
    y = m.add_var(ub=100, cost=1, name="y")
    m.add_constr(x >= 3)
    m.add_pwl(x, y, [0, 2, 4, 6], [0, 4, 4.5, 10], name="cost")
    return m


def abs_model():
    """max r, r = |x|, -5 <= x <= 4: 5 (needs a binary: the LP epigraph would be unbounded)."""
    m = rs.Model(sense="max")
    x = m.add_var(lb=-5, ub=4, name="x")
    r = m.add_var(ub=100, cost=1, name="r")
    m.add_abs(r, x, name="absx")
    return m


def minmax_model():
    """min r - s, r = max(a, b, 1), s = min(a, b): r - s = |a - b| >= 0, and s needs a binary."""
    m = rs.Model()
    a = m.add_var(lb=2, ub=5, name="a")
    b = m.add_var(lb=1.5, ub=5, name="b")
    r = m.add_var(ub=10, cost=1, name="r")
    s = m.add_var(ub=10, cost=-1, name="s")
    m.add_max(r, [a, b], constant=1, name="mx")
    m.add_min(s, [a, b], name="mn")
    return m


def logic_model():
    """a = b = 1 forced; z = and(a, b), w = or(a, b), v = not(a); min -z - w - v: -2."""
    m = rs.Model()
    a, b, z, w, v = (m.add_var(binary=True, name=n) for n in "abzwv")
    m.add_constr(a == 1)
    m.add_constr(b == 1)
    m.set_cost([z, w, v], [-1, -1, -1])
    m.add_and(z, [a, b])
    m.add_or(w, [a, b])
    m.add_not(v, a)
    return m


MPS_KINDS = {"indicator", "sos1", "sos2"}  # the MPS sections the writer has

CASES = [
    ("indicator", lambda: indicator_model()[0], 8.0),
    ("sos1", sos1_model, 3.0),
    ("sos2", sos2_model, 6.0),
    ("pwl", pwl_model, 4.25),
    ("abs", abs_model, 5.0),
    ("minmax", minmax_model, 0.0),
    ("logic", logic_model, -2.0),
]


@pytest.mark.parametrize("label,build,expected", CASES, ids=[c[0] for c in CASES])
def test_golden_path(label, build, expected):
    m = build()
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(expected, abs=1e-6)
    report = sol.general_constraints
    assert report["verified"] is True and report["violations"] == []
    assert len(sol.x) == m.num_vars  # the solution is cut back to the user's columns
    assert len(sol.row_activity) == m.num_constrs
    assert m.check(sol)["verdict"] == "feasible"


@pytest.mark.parametrize("label,build,expected", CASES, ids=[c[0] for c in CASES])
@pytest.mark.parametrize("ext", ["lp", "mps"])
def test_round_trip_through_the_cli(tmp_path, label, build, expected, ext):
    """Build in Python, write the file, `rusolve mip` reads it: the same objective."""
    m = build()
    if ext == "mps" and label not in MPS_KINDS:
        with pytest.raises(rs.ModelError, match="write .lp instead"):
            m.write(str(tmp_path / "x.mps"))
        return
    path = tmp_path / f"{label}.{ext}"
    m.write(str(path))
    assert cli_objective(path) == pytest.approx(m.solve().objective, abs=1e-6)


def test_minmax_values():
    sol = minmax_model().solve()
    a, b, r, s = sol.x
    assert r == pytest.approx(max(a, b, 1.0), abs=1e-6)
    assert s == pytest.approx(min(a, b), abs=1e-6)


def test_pwl_route_and_x_domain():
    sol = pwl_model().solve()
    assert sol.general_constraints["routes"][0]["route"] == "incremental"
    assert sol.x[0] == pytest.approx(3.0, abs=1e-6)
    assert sol.x[1] == pytest.approx(4.25, abs=1e-6)


def test_convex_pwl_is_an_lp():
    m = rs.Model()
    x = m.add_var(ub=10)
    y = m.add_var(ub=100, cost=1)
    m.add_constr(x >= 6)
    m.add_pwl(x, y, [0, 4, 10], [0, 4, 16])
    sol = m.solve()
    assert sol.objective == pytest.approx(8.0)
    assert sol.general_constraints["routes"][0]["route"] == "epigraph_lp"
    assert sol.general_constraints["added_columns"] == 0


def test_sos_differs_from_no_set():
    m = rs.Model(sense="max")
    x = m.add_vars(3, ub=1)
    m.maximize(5 * x[0] + x[1] + 5 * x[2])
    assert m.solve().objective == pytest.approx(11.0)


# ---------------------------------------------------------------- indicator specifics


def test_indicator_matches_its_big_m_twin():
    m, _, _ = indicator_model()
    assert m.solve().objective == pytest.approx(bigm_twin().solve().objective, abs=1e-7)


def test_indicator_on_one_and_equality():
    """min -5 z + x, z = 1 -> x == 4: z = 1 gives -1, z = 0 gives 0."""
    m = rs.Model()
    z = m.add_var(binary=True, cost=-5)
    x = m.add_var(ub=10, cost=1)
    m.add_indicator(z, 1, x == 4)
    sol = m.solve()
    assert sol.objective == pytest.approx(-1.0)
    assert sol.x[0] == pytest.approx(1.0) and sol.x[1] == pytest.approx(4.0)


def test_indicator_envelope_reports_the_object_not_rows():
    m, _, name = indicator_model()
    sol = m.solve()
    assert name == "need_x"
    env = json.loads(sol.to_json())
    g = env["general_constraints"]
    assert g["objects"] == {"indicator": 1}
    assert g["added_rows"] >= 1  # the linearization is reported as such
    assert env["input"]["rows"] == 1 and env["input"]["cols"] == 3  # the user's model
    assert len(env["solution"]["columns"]) == 3


def test_trickle_flow_is_judged_on_rounded_binaries():
    """b = 1e-7 lets a big-M row hold x = 1, but the checker rounds b to 0 and `b = 0 -> x >= 5` fails."""
    m = rs.Model()
    b = m.add_var(binary=True, cost=3, name="b")
    x = m.add_var(ub=10, cost=1, name="x")
    m.add_constr(x + b >= 0)
    m.add_indicator(b, 0, x >= 5, name="i")
    report = m.check({"b": 1e-7, "x": 1.0})
    assert report["verdict"] != "feasible"
    assert any(f["kind"] == "indicator" and f["name"] == "i" for f in report["findings"])
    assert m.check({"b": 1.0, "x": 0.0})["verdict"] == "feasible"
    sol = m.solve()
    assert sol.objective == pytest.approx(3.0)
    assert sol.general_constraints["verified"] is True


def test_check_reports_on_the_user_objects():
    m = sos1_model()
    r = m.check([1.0, 1.0, 0.0])
    assert r["verdict"] != "feasible"
    assert r["summary"]["general_constraints"]["count"] == 1


def test_gc_objects_added_after_a_solve():
    m = rs.Model(sense="max")
    x = m.add_vars(2, ub=1)
    m.maximize(x[0] + x[1])
    assert m.solve().objective == pytest.approx(2.0)
    m.add_sos1(x)
    assert m.solve().objective == pytest.approx(1.0)
    assert m.num_general_constraints == 1


def test_plain_model_has_no_report():
    m = rs.Model(sense="max")
    x = m.add_var(ub=1, binary=True)
    m.maximize(x)
    sol = m.solve()
    assert sol.general_constraints is None
    assert "general_constraints" not in json.loads(sol.to_json())


# ---------------------------------------------------------------- errors


def code(excinfo):
    return excinfo.value.code


def test_indicator_needs_a_binary():
    m = rs.Model()
    x = m.add_var(ub=10)
    with pytest.raises(rs.ModelError) as e:
        m.add_indicator(x, 1, x >= 2)
    assert code(e) == "RS4009" and "not binary" in str(e.value)
    assert m.num_general_constraints == 0


def test_indicator_value_and_row_form():
    m = rs.Model()
    z = m.add_var(binary=True)
    x = m.add_var(ub=10)
    with pytest.raises(rs.ModelError) as e:
        m.add_indicator(z, 2, x >= 2)
    assert code(e) == "RS4009"
    with pytest.raises(TypeError):
        m.add_indicator(z, 1, 5)


def test_indicator_range_row_is_refused():
    m = rs.Model()
    z = m.add_var(binary=True)
    x = m.add_var(ub=10)
    with pytest.raises(rs.ModelError) as e:
        m.add_indicator(z, 1, x.between(1, 5))
    assert code(e) == "RS4009" and "one-sided" in str(e.value)


def test_unknown_variable():
    m = rs.Model()
    x = m.add_var(ub=10)
    with pytest.raises(rs.ModelIndexError) as e:
        m.add_pwl(x, 99, [0, 1], [0, 1])
    assert code(e) == "RS4002"
    other = rs.Model()
    y = other.add_var(binary=True)
    with pytest.raises(rs.ModelError) as e:
        m.add_indicator(y, 1, x >= 1)
    assert code(e) == "RS4003"


def test_pwl_breakpoint_errors():
    m = rs.Model()
    x = m.add_var(ub=10)
    y = m.add_var(ub=10)
    for xs, ys, text in [
        ([0, 2, 1], [0, 1, 2], "sorted"),
        ([0], [0], "at least two"),
        ([0, 1, 2], [0, 1], "one y per breakpoint"),
        ([0, float("nan")], [0, 1], "finite"),
    ]:
        with pytest.raises(rs.ModelError, match=text) as e:
            m.add_pwl(x, y, xs, ys)
        assert code(e) == "RS4009"
    assert m.num_general_constraints == 0


def test_pwl_jump_is_a_discontinuity():
    """f jumps from 1 to 3 at x = 2: at x = 2 either value is allowed."""
    for cost, expected in ((1, 1.0), (-1, -3.0)):
        m = rs.Model()
        x = m.add_var(ub=10)
        y = m.add_var(ub=10, cost=cost)
        m.add_constr(x >= 2 if cost > 0 else x <= 2)
        m.add_pwl(x, y, [0, 2, 2, 4], [0, 1, 3, 4])
        assert m.solve().objective == pytest.approx(expected)


def test_sos_errors():
    m = rs.Model()
    x = m.add_vars(3, ub=1)
    with pytest.raises(rs.ModelError, match="3 variables but 2 weights") as e:
        m.add_sos1(x, weights=[1, 2])
    assert code(e) == "RS4009"
    with pytest.raises(rs.ModelError, match="must differ"):
        m.add_sos2(x, weights=[1, 1, 2])
    with pytest.raises(rs.ModelError, match="appears twice"):
        m.add_sos1([x[0], x[0]])


def test_function_errors():
    m = rs.Model()
    a, b = m.add_vars(2, binary=True)
    c = m.add_var(ub=3)
    with pytest.raises(rs.ModelError, match="not binary") as e:
        m.add_and(c, [a, b])
    assert code(e) == "RS4009"
    with pytest.raises(rs.ModelError, match="exactly one"):
        m.add_abs(c, [a, b])
    with pytest.raises(rs.ModelError, match="at least one"):
        m.add_max(c, [])


def test_duplicate_name():
    m = rs.Model()
    a, b = m.add_vars(2, binary=True)
    m.add_or(b, [a], name="f")
    with pytest.raises(rs.ModelError, match="already used") as e:
        m.add_and(b, [a], name="f")
    assert code(e) == "RS4009"


def test_bounds_changed_after_adding_are_rechecked_at_solve():
    m, (z, x, y), _ = indicator_model()
    m.set_integer(z, False)
    with pytest.raises(rs.ModelError) as e:
        m.solve()
    assert code(e) == "RS4009" and "not binary" in str(e.value)


def test_diagnostics_that_would_ignore_the_objects_are_refused():
    m, _, _ = indicator_model()
    for call in (m.explain, m.iis, m.diagnose, m.ranging, m.feasrelax):
        with pytest.raises(rs.ModelError) as e:
            call()
        assert code(e) == "RS4010"
    with pytest.raises(rs.ModelError) as e:
        m.solve(callback=object())
    assert code(e) == "RS4010"
