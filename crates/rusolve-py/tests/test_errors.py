"""Typed exceptions with `.diagnostics`, and the optional capture of the Python source line.

The diagnostics are the `rusolve.diagnostics/1` documents of docs/design/diagnostics.md: each
exception carries the registered `code`, and `location` is the line of the caller's Python that
raised (a `kind: "python"` location). Run after `maturin develop --release`.
"""
import doctest
import json
import pathlib
import pickle
import time

import pytest

import rusolve as rs
import rusolve.errors as errors

ROOT = pathlib.Path(__file__).resolve().parents[3]


def infeasible_model(capture=False):
    m = rs.Model(name="bad")
    m.capture_source = capture
    x = m.add_var(name="x")
    m.add_constr(x >= 3, name="lo")
    m.add_constr(x <= 1, name="hi")
    m.minimize(x)
    return m


def test_hierarchy_keeps_the_builtin_bases():
    assert issubclass(rs.ModelError, (rs.RusolveError, ValueError))
    assert issubclass(rs.ModelIndexError, (rs.ModelError, IndexError))
    assert issubclass(rs.InfeasibleError, rs.SolveError)
    assert issubclass(rs.UnboundedError, rs.SolveError)
    assert issubclass(rs.LimitError, rs.SolveError)
    assert issubclass(rs.NoSolutionError, (rs.SolveError, ValueError))
    assert rs.errors is errors


def test_model_error_has_code_message_help_and_python_location():
    m = rs.Model()
    with pytest.raises(rs.ModelError) as info:
        m.add_var(lb=float("nan"), name="x")  # the line the location must point at
    e = info.value
    assert isinstance(e, ValueError)
    assert e.code == "RS4001"
    d = e.diagnostics[0]
    assert d["severity"] == "error" and d["message"] and d["help"]
    loc = e.location
    assert loc["kind"] == "python"
    assert loc["path"].endswith("test_errors.py")
    assert loc["source_line"] == 'm.add_var(lb=float("nan"), name="x")  # the line the location must point at'
    assert loc["line"] > 1


def test_mixed_models_is_rs4003():
    a, b = rs.Model(), rs.Model()
    x = a.add_var()
    with pytest.raises(rs.ModelError) as info:
        b.add_constr(x <= 1)
    assert info.value.code == "RS4003"
    assert info.value.location["source_line"] == "b.add_constr(x <= 1)"


def test_index_error_is_rs4002_and_an_indexerror():
    m = rs.Model()
    m.add_vars(2)
    with pytest.raises(IndexError) as info:
        m.add_constr({5: 1.0} <= 1) if False else m.add_constr({5: 1.0}, "<=", 1)
    e = info.value
    assert isinstance(e, rs.ModelIndexError)
    assert e.code == "RS4002"
    assert e.location["kind"] == "python"


def test_named_constraint_keeps_the_type_and_the_code():
    m = rs.Model()
    m.add_vars(2)
    with pytest.raises(rs.ModelIndexError) as info:
        m.add_constr({9: 1.0}, "<=", 1, name="cap")
    assert "cap" in str(info.value)
    assert info.value.code == "RS4002"


def test_infeasible_raise_for_status():
    sol = infeasible_model().solve()
    with pytest.raises(rs.InfeasibleError) as info:
        sol.raise_for_status()  # location line
    e = info.value
    assert e.code in {"RS4004", "RS4005"}  # presolve may report InfeasibleOrUnbounded
    assert e.location["source_line"] == "sol.raise_for_status()  # location line"
    assert any("explain" in h for h in e.help)


def test_infeasible_status_is_rs4004_without_presolve():
    sol = infeasible_model().solve(presolve=False)
    assert sol.status == "Infeasible"
    with pytest.raises(rs.InfeasibleError) as info:
        sol.raise_for_status()
    assert info.value.code == "RS4004"


def test_unbounded_raise_for_status():
    m = rs.Model()
    x = m.add_var(name="x")
    m.maximize(x)
    sol = m.solve(presolve=False)
    assert sol.status == "Unbounded"
    with pytest.raises(rs.UnboundedError) as info:
        sol.raise_for_status()
    assert info.value.code == "RS4005"


def test_limit_error_is_rs4006():
    d = rs.errors  # the class exists; a limit status needs a real limit, so build it from a solve
    m = rs.Model(sense="max")
    x = m.add_vars(40, binary=True)
    import random

    rnd = random.Random(1)
    for _ in range(30):
        m.add_constr(sum(rnd.randint(1, 9) * v for v in x) <= 90)
    m.maximize(sum(rnd.randint(1, 9) * v for v in x))
    sol = m.solve(time_limit=1e-4, presolve=False)
    if sol.status == "Optimal":
        pytest.skip("solved inside the time limit")
    with pytest.raises(d.LimitError) as info:
        sol.raise_for_status()
    assert info.value.code == "RS4006"


def test_optimal_returns_the_solution():
    m = rs.Model(sense="max")
    x = m.add_var(ub=4)
    m.maximize(x)
    sol = m.solve()
    assert sol.raise_for_status() is sol


def test_reading_values_without_a_solution_is_rs4007():
    sol = infeasible_model().solve(presolve=False)
    with pytest.raises(rs.NoSolutionError) as info:
        sol.objective_value if False else sol[0]
    e = info.value
    assert isinstance(e, ValueError)
    assert e.code == "RS4007"
    assert e.location["kind"] == "python"


def test_diagnostics_follow_the_schema():
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((ROOT / "docs/schema/rusolve-diagnostics-1.schema.json").read_text())
    m = rs.Model()
    with pytest.raises(rs.ModelError) as info:
        m.add_var(lb=float("nan"))
    doc = {"schema_version": "rusolve.diagnostics/1", "diagnostics": info.value.diagnostics, "exit_code": 1}
    doc = {k: v for k, v in doc.items() if k in schema.get("properties", {})}
    jsonschema.validate(doc, schema)


def test_registered_codes_explain():
    import subprocess, shutil

    exe = pathlib.Path(rs.__file__).parent / "_bin" / "rusolve-cli"
    exe = exe if exe.exists() else ROOT / "target/release/rusolve-cli"
    if not exe.exists():
        pytest.skip("no CLI binary")
    for code in ["RS4001", "RS4002", "RS4003", "RS4004", "RS4005", "RS4006", "RS4007"]:
        out = subprocess.run([str(exe), "--explain", code], capture_output=True, text=True)
        assert out.returncode == 0 and code in out.stdout


def test_pickle_and_copy_keep_diagnostics():
    m = rs.Model()
    with pytest.raises(rs.ModelError) as info:
        m.add_var(lb=float("nan"))
    e2 = pickle.loads(pickle.dumps(info.value))
    assert type(e2) is rs.ModelError and e2.code == "RS4001" and e2.diagnostics == info.value.diagnostics


def test_errors_module_doctest():
    r = doctest.testmod(errors)
    assert r.failed == 0 and r.attempted > 0


# ------------------------------------------------------------------ source capture


def test_capture_is_off_by_default():
    m = rs.Model()
    assert m.capture_source is False
    x = m.add_var()
    r = m.add_constr(x <= 1)
    assert m.source_of(r) is None


def test_source_of_names_the_adding_line():
    m = rs.Model()
    m.capture_source = True
    x = m.add_vars(3)
    r0 = m.add_constr(x[0] <= 1)
    rs_ = m.add_constrs([x[1] <= 1, x[2] <= 2])
    s0 = m.source_of(r0)
    assert s0["source_line"] == "r0 = m.add_constr(x[0] <= 1)"
    assert s0["path"].endswith("test_errors.py")
    s1 = m.source_of(rs_[1])
    assert s1["source_line"] == "rs_ = m.add_constrs([x[1] <= 1, x[2] <= 2])"
    assert s1["line"] == m.source_of(rs_[0])["line"]
    assert m.source_of(99) is None


def test_explain_puts_the_source_on_each_row():
    m = infeasible_model(capture=True)
    e = m.explain()
    rows = e["explanation"]["iis"]["rows"]
    assert {r["name"] for r in rows} == {"lo", "hi"}
    lines = {r["name"]: r["source"]["source_line"] for r in rows}
    assert lines["lo"] == 'm.add_constr(x >= 3, name="lo")'
    assert lines["hi"] == 'm.add_constr(x <= 1, name="hi")'
    assert all(r["source"]["path"].endswith("test_errors.py") for r in rows)


def test_explain_has_no_source_when_capture_is_off():
    rows = infeasible_model().explain()["explanation"]["iis"]["rows"]
    assert all("source" not in r for r in rows)


def test_capture_cost_is_measured_and_bounded():
    def build(capture, n=20000):
        m = rs.Model()
        m.capture_source = capture
        x = m.add_vars(2)
        t = time.perf_counter()
        for _ in range(n):
            m.add_constr(x[0] + x[1] <= 1)
        return (time.perf_counter() - t) / n

    off = min(build(False) for _ in range(3))
    on = min(build(True) for _ in range(3))
    print(f"add_constr: {off * 1e6:.2f} us off, {on * 1e6:.2f} us on")
    assert on < 50e-6  # opt-in; the measured cost is in the README


def test_infeasible_error_carries_the_cli_iis_diagnostic():
    # The same RS2010 the CLI reports, with the conflicting rows, named Python calls as help
    m = rs.Model()
    x, y = m.add_var(name="x"), m.add_var(name="y")
    m.add_constr(x + y >= 3, name="need3")
    m.add_constr(x <= 1, name="capx")
    m.add_constr(y <= 1, name="capy")
    m.minimize(x + y)
    with pytest.raises(rs.InfeasibleError) as info:
        m.solve().raise_for_status()
    e = info.value
    assert e.code == "RS4004"
    d = next(d for d in e.diagnostics if d["code"] == "RS2010")
    text = json.dumps(d)
    assert all(n in text for n in ("need3", "capx", "capy"))
    helps = " ".join(d.get("help", []))
    assert "Model.iis()" in helps and "Model.feasrelax()" in helps and "rusolve explain" not in helps


def test_unbounded_error_carries_the_ray():
    m = rs.Model()
    x, y = m.add_var(name="x"), m.add_var(name="y")
    m.add_constr(x - y <= 1, name="r")
    m.maximize(x + y)
    with pytest.raises(rs.UnboundedError) as info:
        m.solve().raise_for_status()
    e = info.value
    assert e.code in ("RS4005",)
    d = next(d for d in e.diagnostics if d["code"] == "RS2011")
    text = json.dumps(d)
    assert "ray" in text and "x" in text and "y" in text
    assert "rusolve explain" not in text and "rusolve lint" not in text


def test_a_raising_lazy_callback_reports_interrupted():
    import sys

    sys.path.insert(0, str(pathlib.Path(__file__).parent))
    import hooks_common as hc

    class Boom:
        def lazy_constraints(self, ctx):
            raise RuntimeError("x")

    m = hc.tsp_model(hc.cities(9, 4))
    with pytest.raises(rs.CallbackError) as info:
        m.solve(callback=Boom(), threads=1)
    assert info.value.solution.status == "Interrupted"


def test_qp_unbounded_has_no_hessian_blind_ray():
    # min -x - y + 0.5*x^2 on free x, y: unbounded only along y; an LP-only ray could be along x
    m = rs.Model()
    x, y = m.add_var(lb=-1e30, name="x"), m.add_var(lb=-1e30, name="y")
    m.minimize(-x - y)
    m.set_quadratic_objective([[1.0, 0.0], [0.0, 0.0]])
    sol = m.solve()
    assert sol.status in ("Unbounded", "InfeasibleOrUnbounded"), sol.status
    with pytest.raises(rs.UnboundedError) as info:
        sol.raise_for_status()
    assert all(d["code"] != "RS2011" for d in info.value.diagnostics)


def test_qp_infeasible_still_carries_rs2010():
    m = rs.Model()
    x = m.add_var(name="x")
    m.add_constr(x >= 3, name="lo")
    m.add_constr(x <= 1, name="hi")
    m.minimize(x)
    m.set_quadratic_objective([[1.0]])
    with pytest.raises(rs.InfeasibleError) as info:
        m.solve().raise_for_status()
    d = next(d for d in info.value.diagnostics if d["code"] == "RS2010")
    assert "lo" in json.dumps(d) and "hi" in json.dumps(d)
