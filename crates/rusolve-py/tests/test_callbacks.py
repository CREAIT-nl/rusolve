"""Python callbacks (`Model.solve(callback=...)`) and `Model.diagnose()`.

* lazy-constraint TSP in Python equals brute force on the 8 instances of
  `crates/rusolve-ekk/tests/lazy_cuts.rs`;
* a Python cover separator, the same cuts as model rows and the plain solve agree (180);
* `ctx.interrupt()` stops the solve; a hook that raises gives `rusolve.CallbackError` (RS4008);
* Python branching-score and node-selection prototypes change the tree, not the optimum;
* `Model.diagnose()` equals `rusolve diagnose --format json`.
Run after `maturin develop --release`: `pytest crates/rusolve-py/tests`.
"""
import json
import pathlib
import subprocess
import sys

import pytest

import rusolve as rs

ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "crates/rusolve-py/examples"))
import hooks_common as hc  # noqa: E402

TSP_CASES = [(6, 1), (7, 2), (8, 3), (9, 4), (9, 5), (9, 6), (10, 7), (10, 8)]


@pytest.mark.parametrize("n,seed", TSP_CASES)
def test_lazy_tsp_matches_brute_force(n, seed):
    pts = hc.cities(n, seed)
    brute = hc.brute_force_tour(pts)
    m = hc.tsp_model(pts)
    rows = m.num_constrs
    sol = m.solve(callback=hc.SubtourElimination(n), threads=1, options={"mip_rel_gap": 0.0, "mip_abs_gap": 0.0})
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(brute, abs=1e-6)
    assert len(hc.components(list(sol.x), n)) == 1, "the accepted solution must be one tour"
    st = m.callback_stats
    assert st["lazy"]["accepted"] and st["errors"] == 0
    assert st["calls"]["lazy_constraints"] == st["lazy"]["rounds"]
    # the lazy rows lived in the solver for this solve only
    assert m.num_constrs == rows
    assert len(sol.row_activity) == rows


def test_lazy_tsp_needed_rows_and_brute_force_agree_for_small_n():
    pts = hc.cities(7, 2)
    assert hc.brute_force_tour(pts) == pytest.approx(hc.brute_force_permutations(pts), abs=1e-9)
    # without the callback the degree-2 relaxation is shorter or equal (subtours are allowed)
    plain = hc.tsp_model(hc.cities(9, 4)).solve(threads=1).objective
    assert plain <= hc.brute_force_tour(hc.cities(9, 4)) + 1e-9


def test_lazy_model_can_be_solved_again_and_edited():
    pts = hc.cities(9, 4)
    m = hc.tsp_model(pts)
    a = m.solve(callback=hc.SubtourElimination(9), threads=1).objective
    b = m.solve(callback=hc.SubtourElimination(9), threads=1).objective
    assert a == pytest.approx(b, abs=1e-9)
    x = m.var(0)
    m.add_constr(x <= 1)  # a row added after a lazy solve lands after the model's own rows
    assert m.solve(callback=hc.SubtourElimination(9), threads=1).objective == pytest.approx(a, abs=1e-9)


def test_lazy_round_limit_is_reported_not_hidden():
    m = hc.tsp_model(hc.cities(9, 4))
    sol = m.solve(callback=hc.SubtourElimination(9), lazy_rounds=1, threads=1)
    assert sol.status == "LazyRoundLimit" and sol.x is None and sol.objective is None
    assert m.callback_stats["lazy"]["accepted"] is False


# ---------------------------------------------------------------- user cuts


def test_python_cover_separator_equals_cuts_in_the_model_and_the_plain_solve():
    values, weights, cap = hc.knapsack_data()
    covers = hc.minimal_covers(weights, cap)
    assert len(covers) > 20
    plain = hc.knapsack_model().solve(presolve=False, threads=1, options=hc.CUT_OPTIONS).objective
    in_model = hc.knapsack_model(covers).solve(presolve=False, threads=1, options=hc.CUT_OPTIONS).objective
    m = hc.knapsack_model()
    hooked = m.solve(callback=hc.CoverSeparator(covers), threads=1, options=hc.CUT_OPTIONS)
    assert plain == in_model == hooked.objective == 180.0
    st = m.callback_stats
    assert st["calls"]["separate"] > 0 and st["cuts_accepted"] > 0 and st["errors"] == 0
    assert st["cuts_rejected"] == 0


def test_separator_cuts_as_constraints_and_dict_coefficients():
    values, weights, cap = hc.knapsack_data()
    covers = hc.minimal_covers(weights, cap)

    class Mixed:
        def __init__(self, m):
            self.x = [m.var(j) for j in range(14)]

        def separate(self, ctx):
            out = []
            for k, cv in enumerate(covers):
                if k % 3 == 0:
                    out.append(rs.quicksum(self.x[i] for i in cv) <= len(cv) - 1)  # a Constraint
                elif k % 3 == 1:
                    out.append(({i: 1.0 for i in cv}, len(cv) - 1.0))  # (dict, rhs)
                else:
                    out.append(([(i, 1.0) for i in cv], len(cv) - 1.0))  # ([(col, val)], rhs)
            return out

    m = hc.knapsack_model()
    assert m.solve(callback=Mixed(m), threads=1, options=hc.CUT_OPTIONS).objective == 180.0
    assert m.callback_stats["cuts_accepted"] > 0


def test_malformed_cut_is_rejected_by_the_gate_not_the_solve():
    class Junk:
        def separate(self, ctx):
            return [([999], [1.0], 1.0), ([0, 0], [1.0, 1.0], 1.0)]

    m = hc.knapsack_model()
    assert m.solve(callback=Junk(), threads=1, options=hc.CUT_OPTIONS).objective == 180.0
    assert m.callback_stats["cuts_rejected"] > 0 and m.callback_stats["errors"] == 0


# ---------------------------------------------------------------- events, interrupt, errors


def test_interrupt_stops_the_solve():
    class Stop:
        calls = 0

        def on_progress(self, ctx):
            Stop.calls += 1
            ctx.interrupt()

    m = hc.multi_knapsack()
    sol = m.solve(callback=Stop(), presolve=False, threads=1, options=hc.TREE_OPTIONS)
    assert sol.status == "Interrupted", sol.status
    assert Stop.calls >= 1 and m.callback_stats["interrupted"] is True
    full = hc.multi_knapsack().solve(presolve=False, threads=1, options=hc.TREE_OPTIONS)
    assert full.status == "Optimal" and sol.nodes < full.nodes
    with pytest.raises(rs.LimitError):
        sol.raise_for_status()


def test_interrupt_from_an_improving_solution_event_and_reading_it():
    seen = []

    class Watch:
        def on_improving_solution(self, ctx):
            seen.append((ctx.objective, list(ctx.solution)))
            if len(seen) == 2:
                ctx.interrupt()

    m = hc.multi_knapsack()
    sol = m.solve(callback=Watch(), presolve=False, threads=1, options=hc.TREE_OPTIONS)
    assert len(seen) >= 1
    assert all(len(x) == 28 for _, x in seen)
    if len(seen) >= 2:
        assert sol.status == "Interrupted"


def test_callback_exception_is_a_typed_error_with_a_registered_code():
    class Boom:
        def on_progress(self, ctx):
            raise ValueError("bad hook")

    m = hc.multi_knapsack()
    with pytest.raises(rs.CallbackError) as info:
        m.solve(callback=Boom(), presolve=False, threads=1, options=hc.TREE_OPTIONS)
    e = info.value
    assert e.code == "RS4008" and isinstance(e, rs.RusolveError)
    assert isinstance(e.__cause__, ValueError) and "bad hook" in str(e)
    assert "on_progress" in str(e) and e.location["kind"] == "python"
    assert e.solution is not None and e.solution.status == "Interrupted"
    assert m.callback_stats["errors"] == 1
    # the model is usable afterwards, with no callback
    assert m.solve(presolve=False, threads=1, options=hc.TREE_OPTIONS).status == "Optimal"


@pytest.mark.parametrize(
    "hook,body",
    [
        ("separate", "return 7"),  # not a cut list
        ("branching_score", "return [1.0]"),  # wrong length
        ("select_node", "return 10**6"),  # out of range
        ("lazy_constraints", "raise RuntimeError('x')"),
        ("user_solution", "return [0.0]"),  # wrong length
        ("on_improving_solution", "1 / 0"),
    ],
)
def test_every_hook_failure_is_caught(hook, body):
    ns = {}
    exec(f"class C:\n    def {hook}(self, ctx):\n        {body}\n", ns)
    m = hc.multi_knapsack(34, 5)  # a tree that pops nodes from the queue
    with pytest.raises(rs.CallbackError) as info:
        m.solve(callback=ns["C"](), threads=1, options=hc.TREE_OPTIONS)
    assert info.value.code == "RS4008" and hook in str(info.value)


def test_pickle_keeps_the_code():
    import pickle

    e = rs.CallbackError("x", [{"code": "RS4008"}])
    assert pickle.loads(pickle.dumps(e)).code == "RS4008"


def test_a_callback_without_hooks_is_an_error():
    m = hc.knapsack_model()
    with pytest.raises(rs.ModelError):
        m.solve(callback=object())


def test_user_solution_is_checked_and_used():
    m0 = hc.multi_knapsack()
    opt = m0.solve(presolve=False, threads=1, options=hc.TREE_OPTIONS)
    base_nodes = opt.nodes
    x = list(opt.x)

    class Give:
        def user_solution(self, ctx):
            return x

    m = hc.multi_knapsack()
    s = m.solve(callback=Give(), presolve=False, threads=1, options=hc.TREE_OPTIONS)
    assert s.objective == pytest.approx(opt.objective)
    assert s.nodes < base_nodes
    assert m.callback_stats["user_solutions_accepted"] >= 1

    class GiveBad:
        def user_solution(self, ctx):
            return [1.0] * 28  # violates every row

    m = hc.multi_knapsack()
    s = m.solve(callback=GiveBad(), presolve=False, threads=1, options=hc.TREE_OPTIONS)
    assert s.objective == pytest.approx(opt.objective)
    assert m.callback_stats["user_solutions_rejected"] >= 1


def test_cut_pool_is_readable():
    pools = []

    class Pool:
        def on_cut_pool(self, ctx):
            pools.append(ctx.cut_pool)

    m = hc.knapsack_model()
    m.solve(callback=Pool(), threads=1)
    # the callback may or may not fire on a tiny model; when it does the rows are well formed
    for pool in pools:
        for idx, val, lo, up in pool:
            assert len(idx) == len(val)


# ---------------------------------------------------------------- Python prototypes of steering hooks


def test_python_branching_and_node_selection_change_the_tree_not_the_optimum():
    base = hc.multi_knapsack().solve(presolve=False, threads=1, options=hc.TREE_OPTIONS)
    assert base.status == "Optimal" and base.nodes > 10
    nodes = {}
    for name, hook in [
        ("high", lambda ctx: [float(c) for c, _ in ctx.candidates]),
        ("low", lambda ctx: [-float(c) for c, _ in ctx.candidates]),
    ]:
        class B:
            branching_score = staticmethod(lambda ctx, _h=hook: _h(ctx))

        m = hc.multi_knapsack()
        s = m.solve(callback=B(), threads=1, options=hc.TREE_OPTIONS)
        assert s.status == "Optimal" and s.objective == pytest.approx(base.objective)
        assert m.callback_stats["branching_overrides"] > 0
        nodes[name] = s.nodes
    assert nodes["high"] != nodes["low"] or nodes["high"] != base.nodes

    class Deep:
        def select_node(self, ctx):
            return max(range(len(ctx.open_nodes)), key=lambda k: ctx.open_nodes[k][2])

    # the 28-item tree is solved without queue pops under the default profile; a slightly larger
    # one pops nodes (HighsMipSolver::run installNodes), so select_node is consulted
    plain = hc.multi_knapsack(34, 5).solve(presolve=False, threads=1, options=hc.TREE_OPTIONS)
    m = hc.multi_knapsack(34, 5)
    s = m.solve(callback=Deep(), threads=1, options=hc.TREE_OPTIONS)
    assert s.status == "Optimal" and s.objective == pytest.approx(plain.objective)
    assert m.callback_stats["calls"]["select_node"] > 0


def test_hooks_that_need_model_columns_turn_presolve_off():
    m = hc.knapsack_model()

    class Sep:
        def separate(self, ctx):
            assert len(ctx.solution) == 14  # the model's own columns
            return []

    for kw in ({"presolve": True}, {"presolve": "choose"}, {}):
        s = m.solve(callback=Sep(), threads=1, **kw)
        assert s.status == "Optimal"
        assert dict(s.to_dict()["run"]["options"])["presolve"] == "off"
    # the next solve without a callback neither records nor applies the override
    s = m.solve(threads=1)
    assert s.objective == 180.0
    assert "presolve" not in dict(s.to_dict()["run"]["options"])


def test_no_callback_keeps_the_plain_path():
    m = hc.knapsack_model()
    assert m.solve(threads=1).objective == 180.0
    assert m.callback_stats is None


# ---------------------------------------------------------------- Model.diagnose


def find_cli():
    exe = pathlib.Path(rs.__file__).parent / "_bin" / "rusolve-cli"
    for c in (exe, ROOT / "target/release/rusolve-cli"):
        if c.exists():
            return c
    return None


def infeasible_model():
    m = rs.Model(name="infeasible")
    x = m.add_var(lb=0, ub=10, name="x")
    y = m.add_var(lb=0, ub=10, name="y")
    m.add_constr(x + y >= 30, name="demand")
    m.add_constr(x + y <= 12, name="capacity")
    m.minimize(x + y)
    return m


def test_diagnose_returns_the_cli_object_for_an_infeasible_model():
    d = infeasible_model().diagnose()
    assert d["schema_version"] == "rusolve.diagnose/1"
    assert d["status"] == "Infeasible" and d["verdict"] in {"formulation", "data", "unclear"}
    assert d["kind"] == "lp" and d["scores"].keys() == {"data", "formulation"}
    assert json.dumps(d)  # plain JSON


def test_diagnose_of_a_feasible_model_says_nothing_to_diagnose():
    d = hc.knapsack_model().diagnose()
    assert d["verdict"] == "none" and "nothing to diagnose" in d["summary"]
    assert d["kind"] == "relaxation"


def test_diagnose_equals_the_cli(tmp_path):
    exe = find_cli()
    if exe is None:
        pytest.skip("no CLI binary")
    for build in (infeasible_model, hc.knapsack_model):
        m = build()
        path = tmp_path / f"{build.__name__}.mps"
        m.write(str(path))
        proc = subprocess.run([str(exe), "diagnose", str(path), "--format", "json"], capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        cli = json.loads(proc.stdout)
        py = m.diagnose()
        for d in (cli, py):
            d["input"].pop("path")
        assert py == cli
        assert py["verdict"] == cli["verdict"]


def test_diagnose_refuses_a_quadratic_objective():
    m = rs.Model()
    x = m.add_var(name="x")
    m.set_quadratic_objective({(0, 0): 2.0})
    with pytest.raises(rs.ModelError) as info:
        m.diagnose()
    assert info.value.code == "RS4001"
