"""Acceptance tests for the rusolve Python bindings.

Run after `maturin develop --release` (see ../README.md): `pytest crates/rusolve-py/tests`.
"""
import math
import os
import pathlib
import subprocess

import pytest

import rusolve

ROOT = pathlib.Path(__file__).resolve().parents[3]
KNAP = ROOT / "crates/rusolve-ekk/tests/fixtures/s5/e2e_knap.mps"
CLI = pathlib.Path(os.environ.get("CARGO_TARGET_DIR", ROOT / "target")) / "release/rusolve-cli"


def quiet(h):
    # the options `rusolve mip` sets (mip_highs.rs::trace_options)
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("random_seed", 0)
    return h


def example_lp():
    """highspy-style example: min x + y s.t. x + 2y >= 2, 3x + y >= 3, x, y >= 0."""
    inf = rusolve.kHighsInf
    lp = rusolve.HighsLp()
    lp.num_col_ = 2
    lp.num_row_ = 2
    lp.col_cost_ = [1.0, 1.0]
    lp.col_lower_ = [0.0, 0.0]
    lp.col_upper_ = [inf, inf]
    lp.row_lower_ = [2.0, 3.0]
    lp.row_upper_ = [inf, inf]
    lp.a_matrix_.format_ = rusolve.MatrixFormat.kColwise
    lp.a_matrix_.start_ = [0, 2, 4]
    lp.a_matrix_.index_ = [0, 1, 0, 1]
    lp.a_matrix_.value_ = [1.0, 3.0, 2.0, 1.0]
    return lp


def test_highspy_style_lp():
    h = quiet(rusolve.Highs())
    assert h.passModel(example_lp()) == rusolve.HighsStatus.kOk
    assert h.run() == rusolve.HighsStatus.kOk
    assert h.getModelStatus() == rusolve.HighsModelStatus.kOptimal
    sol = h.getSolution()
    # optimum at the intersection x + 2y = 2, 3x + y = 3: x = 0.8, y = 0.6
    assert sol.col_value == pytest.approx([0.8, 0.6])
    assert sol.row_value == pytest.approx([2.0, 3.0])
    assert len(sol.col_dual) == 2 and len(sol.row_dual) == 2
    assert h.getInfo().objective_function_value == pytest.approx(1.4)
    assert h.getObjectiveValue() == pytest.approx(1.4)


def test_flat_pass_model_and_modification():
    inf = rusolve.kHighsInf
    h = quiet(rusolve.Highs())
    st = h.passModel(2, 2, 4, rusolve.MatrixFormat.kColwise, rusolve.ObjSense.kMinimize, 0.0,
                     [1.0, 1.0], [0.0, 0.0], [inf, inf], [2.0, 3.0], [inf, inf],
                     [0, 2, 4], [0, 1, 0, 1], [1.0, 3.0, 2.0, 1.0])
    assert st == rusolve.HighsStatus.kOk
    h.run()
    assert h.getObjectiveValue() == pytest.approx(1.4)
    # add a row x + y >= 3 and a column
    assert h.addRow(3.0, inf, 2, [0, 1], [1.0, 1.0]) == rusolve.HighsStatus.kOk
    h.run()
    assert h.getObjectiveValue() == pytest.approx(3.0)
    assert h.addCol(10.0, 0.0, inf, 1, [0], [1.0]) == rusolve.HighsStatus.kOk
    assert h.getNumCol() == 3 and h.getNumRow() == 3
    h.changeColCost(0, 5.0)
    h.run()
    assert h.getModelStatus() == rusolve.HighsModelStatus.kOptimal


def test_options():
    h = rusolve.Highs()
    assert h.setOptionValue("time_limit", 10) == rusolve.HighsStatus.kOk
    assert h.getOptionValue("time_limit") == (rusolve.HighsStatus.kOk, 10.0)
    assert h.setOptionValue("threads", 1) == rusolve.HighsStatus.kOk
    assert h.getOptionValue("threads")[1] == 1
    assert h.setOptionValue("output_flag", False) == rusolve.HighsStatus.kOk
    assert h.getOptionValue("output_flag")[1] is False
    assert h.setOptionValue("presolve", "off") == rusolve.HighsStatus.kOk
    assert h.getOptionValue("presolve")[1] == "off"
    assert h.setOptionValue("no_such_option", 1) == rusolve.HighsStatus.kError
    with pytest.raises(ValueError):
        rusolve.Highs(profile="nope")
    rusolve.Highs(profile="highs-1.15.1")


def cli_csv():
    if not CLI.exists():
        pytest.skip(f"{CLI} not built")
    env = dict(os.environ, RUSOLVE_MIP_HIGHS="1")
    env.pop("RUSOLVE_PROFILE", None)
    out = subprocess.run([str(CLI), "mip", str(KNAP), "--csv"], env=env, check=True,
                         capture_output=True, text=True).stdout.splitlines()
    return dict(zip(out[0].split(","), out[1].split(",")))


def test_read_mip_matches_cli(tmp_path):
    h = quiet(rusolve.Highs())
    assert h.readModel(str(KNAP)) == rusolve.HighsStatus.kOk
    assert h.run() == rusolve.HighsStatus.kOk
    info = h.getInfo()
    assert h.getModelStatus() == rusolve.HighsModelStatus.kOptimal
    assert info.objective_function_value == pytest.approx(-417.0)
    assert info.mip_node_count == 1 and info.simplex_iteration_count == 10
    csv = cli_csv()
    assert csv["status"] == "Optimal"
    assert float(csv["objective"]) == info.objective_function_value
    assert float(csv["bound"]) == info.mip_dual_bound
    assert int(csv["nodes"]) == info.mip_node_count
    assert int(csv["lp_iters"]) == info.simplex_iteration_count
    xs = [float(v) for v in csv["x"].strip("[]").split(";")]
    assert h.getSolution().col_value == xs
    # solution file with the model's names
    path = tmp_path / "knap.sol"
    assert h.writeSolution(str(path)) == rusolve.HighsStatus.kOk
    text = path.read_text()
    assert "Optimal" in text and "x1 1" in text and "k0 " in text


def test_bound_change_and_resolve():
    h = quiet(rusolve.Highs())
    h.readModel(str(KNAP))
    h.run()
    first = h.getObjectiveValue()
    sol = h.getSolution().col_value
    j = sol.index(1.0)
    # forbid the item used in the optimum: the optimum cannot improve
    assert h.changeColBounds(j, 0.0, 0.0) == rusolve.HighsStatus.kOk
    h.run()
    assert h.getModelStatus() == rusolve.HighsModelStatus.kOptimal
    assert h.getObjectiveValue() >= first - 1e-9
    assert h.getSolution().col_value[j] == 0.0
    assert not math.isnan(h.getInfo().mip_dual_bound)
    # restoring the bound recovers the original optimum
    h.changeColBounds(j, 0.0, 1.0)
    h.run()
    assert h.getObjectiveValue() == pytest.approx(first)
