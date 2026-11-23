"""The diagnostics of the CLI as Model methods, and `Solution.to_json()` (the envelope).

Each method returns the object the CLI prints for `--format json`; the tests check the shape
(`schema_version`), the content on a model whose answer is known, and that the call does not
change the model. Run after `maturin develop --release`: `pytest crates/rusolve-py/tests`.
"""
import json
import os
import pathlib
import subprocess

import pytest

import rusolve as rs

ROOT = pathlib.Path(__file__).resolve().parents[3]
SCHEMAS = ROOT / "docs/schema"


def schema(name):
    return json.loads((SCHEMAS / f"rusolve-{name}-1.schema.json").read_text())


def validate(obj, name):
    jsonschema = pytest.importorskip("jsonschema")
    jsonschema.validate(obj, schema(name))


def knapsack():
    m = rs.Model(sense="max", name="knap")
    x = m.add_vars(4, binary=True, names="x")
    m.add_constr(3 * x[0] + 4 * x[1] + 5 * x[2] + 6 * x[3] <= 10, name="weight")
    m.maximize(4 * x[0] + 5 * x[1] + 6 * x[2] + 8 * x[3])
    return m, x


def infeasible():
    m = rs.Model(name="bad")
    x = m.add_var(name="x")
    y = m.add_var(name="y")
    m.add_constr(x + y >= 3, name="lo")
    m.add_constr(x + y <= 1, name="hi")
    m.add_constr(x <= 5, name="unrelated")
    m.minimize(x + y)
    return m


# ---------------------------------------------------------------- Solution.to_json


def test_to_json_is_the_result_envelope():
    m, _ = knapsack()
    sol = m.solve(time_limit=10)
    env = json.loads(sol.to_json())
    validate(env, "result")
    assert env["schema_version"] == "rusolve.result/1"
    assert env["status"] == "Optimal" and sol.status == "Optimal"
    assert env["objective"] == sol.objective == 13.0
    assert env["input"]["path"] == "knap"
    assert (env["input"]["rows"], env["input"]["cols"], env["input"]["integer_cols"]) == (1, 4, 4)
    assert env["class"] == "mip" and env["has_solution"] is True
    assert env["run"]["options"] == {"time_limit": "10"}
    assert [c["name"] for c in env["solution"]["columns"]] == ["x[0]", "x[1]", "x[2]", "x[3]"]
    assert [c["value"] for c in env["solution"]["columns"]] == list(sol.x)


def test_to_json_options_and_dict():
    m, _ = knapsack()
    sol = m.solve()
    assert "solution" not in json.loads(sol.to_json(include_solution=False))
    assert sol.to_dict() == json.loads(sol.to_json())
    assert "\n  " in sol.to_json(indent=2)


def test_to_json_digest_follows_the_model_and_is_stable():
    m, x = knapsack()
    a = m.solve().to_dict(False)["input"]["sha256"]
    assert m.solve().to_dict(False)["input"]["sha256"] == a
    m.set_cost(x[0], 5.0)
    assert m.solve().to_dict(False)["input"]["sha256"] != a


def test_to_json_of_an_lp_and_an_infeasible_model():
    m = rs.Model()
    x = m.add_var(ub=3)
    m.maximize(x)
    env = m.solve().to_dict()
    validate(env, "result")
    assert env["class"] == "lp" and env["stats"]["nodes"] == 0
    bad = infeasible().solve().to_dict()
    validate(bad, "result")
    assert bad["status"] == "Infeasible" and bad["objective"] is None
    assert "solution" not in bad


def test_to_json_is_a_snapshot_across_resolves():
    m, x = knapsack()
    s1 = m.solve()
    m.set_bounds(x[3], 0, 0)
    s2 = m.solve()
    assert s1.to_dict()["objective"] == 13.0
    assert s2.to_dict()["objective"] == s2.objective != 13.0


def test_to_json_names_follow_added_columns():
    m = rs.Model()
    a = m.add_var(ub=1, name="a")
    m.maximize(a)
    m.solve()
    b = m.add_var(ub=1, name="b")
    m.maximize(a + b)
    env = m.solve().to_dict()
    assert [c["name"] for c in env["solution"]["columns"]] == ["a", "b"]


def test_to_json_qp_names_the_class():
    m = rs.Model()
    x = m.add_vars(2, lb=0, ub=10)
    m.add_constr(x[0] + x[1] >= 1)
    m.set_quadratic_objective({(0, 0): 2.0, (1, 1): 2.0})
    env = m.solve().to_dict()
    assert env["class"] == "convex_qp" and env["convexity"] == "convex"


# ---------------------------------------------------------------- check


def test_check_accepts_a_solution_and_flags_a_violation():
    m, x = knapsack()
    sol = m.solve()
    rep = m.check(sol)
    validate(rep, "check")
    assert rep["verdict"] == "feasible" and rep["input" if "input" in rep else "model"]
    bad = m.check([1.0, 1.0, 1.0, 1.0])
    assert bad["verdict"] != "feasible"
    assert any(f.get("name") == "weight" or "weight" in json.dumps(f) for f in bad["findings"])


def test_check_takes_dicts_by_var_and_by_name():
    m, x = knapsack()
    by_var = {x[i]: v for i, v in enumerate([1.0, 0.0, 0.0, 1.0])}
    assert m.check(by_var)["verdict"] == "feasible"
    by_name = {f"x[{i}]": v for i, v in enumerate([1.0, 1.0, 1.0, 1.0])}
    assert m.check(by_name)["verdict"] != "feasible"


def test_check_rejects_a_wrong_shape():
    m, x = knapsack()
    with pytest.raises(ValueError, match="expected 4 values"):
        m.check([1.0])
    with pytest.raises(ValueError, match="misses 3 of 4"):
        m.check({x[0]: 1.0})
    other, _ = knapsack()
    with pytest.raises(ValueError, match="another Model"):
        m.check(other.solve())
    bad = infeasible()
    with pytest.raises(ValueError, match="no primal values"):
        bad.check(bad.solve())


def test_check_refuses_a_quadratic_objective():
    m = rs.Model()
    x = m.add_vars(2)
    m.set_quadratic_objective({(0, 0): 1.0})
    with pytest.raises(ValueError, match="quadratic"):
        m.check([0.0, 0.0])


# ---------------------------------------------------------------- iis / explain


def test_explain_names_the_conflicting_rows():
    m = infeasible()
    e = m.explain()
    validate(e, "explain")
    assert e["schema_version"] == "rusolve.explain/1"
    text = json.dumps(e)
    assert "lo" in text and "hi" in text
    assert "unrelated" not in text


def test_iis_is_explain():
    m = infeasible()
    assert m.iis() == m.explain()


def test_explain_does_not_change_the_model():
    m = infeasible()
    before = (m.num_vars, m.num_constrs, m.num_nonzeros)
    m.explain()
    assert (m.num_vars, m.num_constrs, m.num_nonzeros) == before
    assert m.solve().status == "Infeasible"


def test_explain_of_an_unbounded_model():
    m = rs.Model(sense="max")
    x = m.add_var()
    m.add_constr(x >= 1)
    m.maximize(x)
    e = m.explain()
    validate(e, "explain")


# ---------------------------------------------------------------- lint


def test_lint_flags_a_big_m():
    m = rs.Model()
    x = m.add_var(name="flow")
    y = m.add_var(binary=True, name="open")
    m.add_constr(x <= 1e9 * y, name="link")
    m.minimize(x + y)
    rep = m.lint()
    validate(rep, "lint")
    assert rep["findings"], "a 1e9 coefficient next to a 1 is a finding"
    rules = {f["rule"] for f in rep["findings"]}
    only = m.lint(rules=[sorted(rules)[0]])
    assert {f["rule"] for f in only["findings"]} == {sorted(rules)[0]}


def test_lint_rejects_an_unknown_rule():
    m, _ = knapsack()
    with pytest.raises(ValueError, match="unknown rule"):
        m.lint(rules=["no-such-rule"])
    with pytest.raises(ValueError, match="min-severity"):
        m.lint(min_severity="loud")  # type: ignore[arg-type]


# ---------------------------------------------------------------- ranging


def test_ranging_of_an_lp():
    m = rs.Model(sense="max")
    x = m.add_vars(2, ub=4, names="v")
    m.add_constr(x[0] + x[1] <= 5, name="cap")
    m.maximize(x[0] + 2 * x[1])
    rep = m.ranging()
    validate(rep, "ranging")
    assert rep["ranging_valid"] is True
    assert {c["name"] for c in rep["columns"]} == {"v[0]", "v[1]"}
    only = m.ranging(rows=["cap"], cols=["v[0]"])
    assert [c["name"] for c in only["columns"]] == ["v[0]"]


def test_ranging_refuses_a_quadratic_objective():
    m = rs.Model()
    m.add_vars(1)
    m.set_quadratic_objective({(0, 0): 1.0})
    with pytest.raises(ValueError, match="quadratic"):
        m.ranging()


# ---------------------------------------------------------------- feasrelax


def test_feasrelax_relaxes_the_conflict():
    m = infeasible()
    rep = m.feasrelax()
    validate(rep, "feasrelax")
    assert rep["status"] in ("relaxed", "optimal", "Optimal") or rep["verified"] is not None
    assert rep["relaxed"] is not None


def test_feasrelax_modes_and_arguments():
    m = infeasible()
    rep = m.feasrelax(mode="min-count", time_limit=10)
    assert rep["mode"] == "min-count"
    rep = m.feasrelax(global_rhs_penalty=1.0, global_lower_penalty=2.0)
    assert rep["mode"] == "weighted"
    with pytest.raises(ValueError, match="mode"):
        m.feasrelax(mode="nope")


# ---------------------------------------------------------------- write


CLI = pathlib.Path(os.environ.get("CARGO_TARGET_DIR", ROOT / "target")) / "release/rusolve-cli"


def cli_objective(path):
    """The objective `rusolve-cli mip <file> --format json` reports (reads MPS and LP)."""
    if not CLI.exists():
        pytest.skip("rusolve-cli is not built")
    out = subprocess.run([str(CLI), "mip", str(path), "--format", "json", "--threads", "1"],
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out)["objective"]


@pytest.mark.parametrize("ext", ["mps", "lp"])
def test_write_round_trips_through_the_cli(tmp_path, ext):
    m, _ = knapsack()
    path = tmp_path / f"knap.{ext}"
    m.write(str(path))
    assert path.stat().st_size > 0
    assert cli_objective(path) == pytest.approx(13.0)


@pytest.mark.parametrize("ext", ["mps", "lp"])
def test_write_gzip_is_the_plain_file_compressed(tmp_path, ext):
    import gzip

    m, _ = knapsack()
    m.write(str(tmp_path / f"k.{ext}"))
    m.write(str(tmp_path / f"k.{ext}.gz"))
    assert gzip.decompress((tmp_path / f"k.{ext}.gz").read_bytes()) == (tmp_path / f"k.{ext}").read_bytes()


def test_write_keeps_ranges_and_offsets(tmp_path):
    m = rs.Model()
    x = m.add_vars(2, ub=10, names="v")
    m.add_constr(x[0] + x[1], lo=1, hi=4, name="rng")
    m.minimize(x[0] + x[1] + 7)
    path = tmp_path / "r.mps"
    m.write(str(path))
    assert "RANGES" in path.read_text()
    assert cli_objective(path) == pytest.approx(m.solve().objective) == pytest.approx(8.0)


def test_write_quadratic_objective_as_lp(tmp_path):
    m = rs.Model()
    x = m.add_vars(2, lb=1, ub=10, names="v")
    m.set_quadratic_objective({(0, 0): 2.0, (1, 1): 2.0})
    m.minimize(x[0] + x[1])
    path = tmp_path / "q.lp"
    m.write(str(path))
    assert "[" in path.read_text()  # the LP format's quadratic section
    assert cli_objective(path) == pytest.approx(m.solve().objective)
    with pytest.raises(ValueError, match="quadratic objective cannot be written as MPS"):
        m.write(str(tmp_path / "q.mps"))


def test_write_errors(tmp_path):
    m, _ = knapsack()
    with pytest.raises(ValueError, match="unknown output format"):
        m.write(str(tmp_path / "m.txt"))
    s = rs.Model()
    v = s.add_var(lb=1, ub=5)
    s.set_semi(v)
    with pytest.raises(ValueError, match="semi-continuous"):
        s.write(str(tmp_path / "s.mps"))


def test_explain_and_lint_refuse_unsound_quadratic_verdicts():
    # min -x + x^2 with x free is bounded (optimum -0.25); the LP part alone is unbounded.
    m = rs.Model()
    x = m.add_var(lb=float("-inf"))
    m.minimize(-1.0 * x)
    m.set_quadratic_objective({(0, 0): 2.0})
    with pytest.raises(rs.ModelError) as info:
        m.explain()
    assert info.value.code == "RS4001"
    with pytest.raises(rs.ModelError):
        m.iis()
    with pytest.raises(rs.ModelError):
        m.lint(with_solve=True)
    assert m.lint()["schema_version"] == "rusolve.lint/1"

    # feasrelax on an infeasible QP: free x, min -x + x^2, rows y>=3 and y<=1 (true relaxed
    # optimum -0.25; the Hessian-blind LP verdict was Unbounded).
    q = rs.Model()
    x = q.add_var(lb=float("-inf"))
    y = q.add_var()
    q.minimize(-1.0 * x)
    q.set_quadratic_objective({(0, 0): 2.0})
    q.add_constr(y >= 3)
    q.add_constr(y <= 1)
    with pytest.raises(rs.ModelError) as info:
        q.feasrelax()
    assert info.value.code == "RS4001"
    with pytest.raises(rs.ModelError) as info:
        q.explain()
    assert info.value.code == "RS4001"
