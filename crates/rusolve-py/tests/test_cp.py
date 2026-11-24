"""Scheduling and combinatorial objects through the Python API: intervals, no-overlap, cumulative,
all-different, element and table (`Model.add_*`, docs/design/cp-model.md). Every model is small
enough to enumerate in the test, so the optimum of the solve is compared with a brute-force
optimum that does not use rusolve, and the answer is checked on the objects (`Solution.cp`,
`Model.check`)."""

import itertools
import json

import pytest

import rusolve as rs

# ---------------------------------------------------------------- job shop

# three jobs, each one operation on machine 1 and then one on machine 2: (duration on M1, on M2)
JOBS = [(2, 3), (3, 1), (2, 2)]


def job_shop():
    m = rs.Model(name="jobshop")
    horizon = sum(a + b for a, b in JOBS)
    makespan = m.add_var(ub=horizon, integer=True, cost=1, name="makespan")
    on = ([], [])
    for j, (d1, d2) in enumerate(JOBS):
        s1 = m.add_var(ub=horizon, integer=True, name=f"s1_{j}")
        s2 = m.add_var(ub=horizon, integer=True, name=f"s2_{j}")
        a = m.add_interval(s1, d1, name=f"a{j}")
        b = m.add_interval(s2, d2, name=f"b{j}")
        m.add_constr(b.start >= a.end)
        m.add_constr(makespan >= b.end)
        on[0].append(a)
        on[1].append(b)
    m.add_no_overlap(on[0], name="machine1")
    m.add_no_overlap(on[1], name="machine2")
    return m, makespan


def brute_job_shop():
    best = None
    for p1 in itertools.permutations(range(3)):
        for p2 in itertools.permutations(range(3)):
            t, end1 = 0, {}
            for j in p1:
                t += JOBS[j][0]
                end1[j] = t
            t = 0
            for j in p2:
                t = max(t, end1[j]) + JOBS[j][1]
            best = t if best is None else min(best, t)
    return best


def test_job_shop_matches_brute_force_and_is_checked():
    m, makespan = job_shop()
    sol = m.solve()
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(brute_job_shop())
    assert sol.cp["verified"] is True and sol.cp["violations"] == []
    assert sol.cp["engine"] == "cp-linearization"
    assert len(sol.x) == m.num_vars
    assert m.check(sol)["verdict"] == "feasible"
    env = json.loads(sol.to_json())
    assert env["class"] == "cp" and env["engine"] == "cp-linearization"
    assert env["cp"]["verified"] is True


def test_check_reports_an_overlap_on_the_object():
    m, makespan = job_shop()
    sol = m.solve()
    x = list(sol.x)
    # put job 0 and job 1 on machine 1 at the same start: the no_overlap must fail
    names = {m.var_name(m.var(j)): j for j in range(m.num_vars)}
    x[names["s1_1"]] = x[names["s1_0"]]
    x[names["a1_end"]] = x[names["s1_1"]] + JOBS[1][0]
    r = m.check(x)
    assert r["verdict"] != "feasible"
    assert any(f["kind"] == "no_overlap" and f["name"] == "machine1" for f in r["findings"])


# ---------------------------------------------------------------- cumulative (RCPSP)

TASKS = [(2, 1), (1, 2), (2, 1), (1, 1)]  # (duration, demand), capacity 2; task 0 before task 2
CAP = 2
H = 6


def rcpsp():
    m = rs.Model(name="rcpsp")
    mk = m.add_var(ub=H, integer=True, cost=1, name="mk")
    ivs = []
    for i, (d, _) in enumerate(TASKS):
        s = m.add_var(ub=H, integer=True, name=f"s{i}")
        ivs.append(m.add_interval(s, d, name=f"t{i}"))
        m.add_constr(mk >= ivs[-1].end)
    m.add_constr(ivs[2].start >= ivs[0].end)
    m.add_cumulative(ivs, [dem for _, dem in TASKS], CAP, name="crew")
    return m


def brute_rcpsp():
    best = None
    for starts in itertools.product(range(H + 1), repeat=len(TASKS)):
        ends = [s + d for s, (d, _) in zip(starts, TASKS)]
        if max(ends) > H or starts[2] < ends[0]:
            continue
        ok = all(
            sum(dem for s, e, (_, dem) in zip(starts, ends, TASKS) if s <= t < e) <= CAP
            for t in range(H + 1)
        )
        if ok:
            best = max(ends) if best is None else min(best, max(ends))
    return best


def test_rcpsp_matches_brute_force():
    m = rcpsp()
    sol = m.solve()
    assert sol.objective == pytest.approx(brute_rcpsp())
    assert sol.cp["verified"] is True
    assert m.check(sol)["verdict"] == "feasible"
    routes = {r["route"] for r in sol.cp["routes"]}
    assert "cumulative_time_indexed" in routes


def test_cumulative_over_the_horizon_limit_is_refused_with_a_typed_error():
    m = rs.Model()
    s = m.add_var(ub=5000, integer=True, name="s")
    t = m.add_var(ub=5000, integer=True, name="t")
    a, b = m.add_interval(s, 3), m.add_interval(t, 4)
    m.add_cumulative([a, b], 1, 1)
    with pytest.raises(rs.ModelError) as e:
        m.solve()
    assert e.value.code == "RS4013" and "horizon" in str(e.value)


# ---------------------------------------------------------------- optional intervals


def test_optional_interval_is_free_when_absent():
    m = rs.Model(sense="max")
    s = m.add_var(ub=3, integer=True, name="s")
    t = m.add_var(ub=3, integer=True, name="t")
    p = m.add_var(binary=True, cost=5, name="p")
    q = m.add_var(binary=True, cost=4, name="q")
    a = m.add_interval(s, 2, presence=p, name="a")
    b = m.add_interval(t, 2, presence=q, name="b")
    m.add_no_overlap([a, b])
    m.add_constr(s == 0)
    m.add_constr(t == 1)
    sol = m.solve()
    # both starting inside each other's window: only the more valuable one can be present
    assert sol.objective == pytest.approx(5.0)
    assert sol.cp["verified"] is True


# ---------------------------------------------------------------- all-different, element, table


def test_all_different_matches_brute_force():
    cost = [3, 1, 2]
    m = rs.Model()
    x = m.add_vars(3, lb=1, ub=3, integer=True, cost=cost, names="x")
    m.add_all_different(x)
    m.add_constr(x[0] >= 2)
    sol = m.solve()
    want = min(
        sum(c * v for c, v in zip(cost, p))
        for p in itertools.permutations([1, 2, 3])
        if p[0] >= 2
    )
    assert sol.objective == pytest.approx(want)
    assert sol.cp["verified"] is True
    assert m.check([1.0, 1.0, 2.0])["verdict"] != "feasible"


def test_element_matches_brute_force():
    values = [5, 3, 8, 2]
    m = rs.Model(sense="max")
    i = m.add_var(ub=3, integer=True, cost=-1, name="i")
    t = m.add_var(ub=10, integer=True, cost=2, name="t")
    m.add_element(i, values, t)
    m.add_constr(i <= 2)
    sol = m.solve()
    want = max(2 * values[k] - k for k in range(3))
    assert sol.objective == pytest.approx(want)
    assert sol.cp["verified"] is True
    assert m.check({"i": 3, "t": 5})["verdict"] != "feasible"


TUPLES = [(0, 1), (1, 3), (2, 2), (3, 0)]


@pytest.mark.parametrize("negated", [False, True])
def test_table_matches_brute_force(negated):
    m = rs.Model(sense="max")
    x = m.add_var(ub=3, integer=True, cost=2, name="x")
    y = m.add_var(ub=3, integer=True, cost=1, name="y")
    m.add_table([x, y], TUPLES, negated=negated)
    sol = m.solve()
    pts = [(a, b) for a in range(4) for b in range(4) if ((a, b) in TUPLES) != negated]
    want = max(2 * a + b for a, b in pts)
    assert sol.objective == pytest.approx(want)
    assert sol.cp["verified"] is True


def test_infeasible_cp_model_has_no_solution():
    m = rs.Model()
    x = m.add_vars(3, lb=0, ub=1, integer=True, names="x")
    m.add_all_different(x)
    sol = m.solve()
    assert not sol.ok and sol.status == "Infeasible"


# ---------------------------------------------------------------- plain models, errors


def test_plain_model_has_no_cp_report():
    m = rs.Model(sense="max")
    x = m.add_var(ub=1, binary=True)
    m.maximize(x)
    sol = m.solve()
    assert sol.cp is None
    assert "cp" not in json.loads(sol.to_json())


def test_invalid_objects_are_typed_and_not_added():
    m = rs.Model()
    x = m.add_var(ub=5)  # continuous
    with pytest.raises(rs.ModelError) as e:
        m.add_all_different([x])
    assert e.value.code == "RS4011" and "not integer" in str(e.value)
    assert m.num_cp_objects == 0
    s = m.add_var(ub=5, integer=True)
    with pytest.raises(rs.ModelError) as e:
        m.add_interval(s, -1)
    assert e.value.code == "RS4011"
    with pytest.raises(rs.ModelError) as e:
        m.add_table([s], [(1, 2)])
    assert e.value.code == "RS4011" and "tuple" in str(e.value)
    with pytest.raises(rs.ModelError) as e:
        m.add_no_overlap([s])
    assert e.value.code == "RS4011"
    assert m.num_cp_objects == 0


def test_calls_that_would_ignore_the_objects_refuse():
    m, _ = job_shop()
    for call in (m.explain, m.iis, m.diagnose, m.ranging, m.feasrelax):
        with pytest.raises(rs.ModelError) as e:
            call()
        assert e.value.code == "RS4012"
    with pytest.raises(rs.ModelError) as e:
        m.write("/tmp/never.lp")
    assert e.value.code == "RS4012"
    with pytest.raises(rs.ModelError) as e:
        m.solve(callback=lambda ctx: None)
    assert e.value.code == "RS4012"


def test_intervals_of_another_model_are_refused():
    a, b = rs.Model(), rs.Model()
    iv = a.add_interval(a.add_var(ub=5, integer=True), 2)
    with pytest.raises(rs.ModelError) as e:
        b.add_no_overlap([iv])
    assert e.value.code == "RS4003"


# ---------------------------------------------------------------- CP-SAT (mip_engine="cpsat")

# the engine loads the OR-Tools library at run time: these tests skip when it is not installed
needs_cpsat = pytest.mark.skipif(
    not rs.cpsat_info()["available"],
    reason="OR-Tools library not found (pip install 'rusolve[cpsat]')",
)


@needs_cpsat
def test_job_shop_with_cpsat_is_checked_on_the_objects():
    m, _ = job_shop()
    sol = m.solve(options={"mip_engine": "cpsat", "mip_rel_gap": 0.0})
    assert sol.status == "Optimal"
    assert sol.objective == pytest.approx(brute_job_shop())
    assert sol.cp["engine"] == "cp-sat"
    assert sol.cp["verified"] is True and sol.cp["violations"] == []
    assert sol.cp["cpsat"]["status"] == "OPTIMAL"
    assert sol.cp["cpsat"]["ortools_version"].startswith("9.")
    assert len(sol.x) == m.num_vars
    assert m.check(sol)["verdict"] == "feasible"
    env = json.loads(sol.to_json())
    assert env["class"] == "cp" and env["engine"] == "cp-sat"


@needs_cpsat
def test_cpsat_takes_what_the_linearization_refuses():
    m = rs.Model()
    s = m.add_var(ub=5000, integer=True, name="s")
    t = m.add_var(ub=5000, integer=True, name="t")
    a, b = m.add_interval(s, 3), m.add_interval(t, 4)
    m.add_cumulative([a, b], 1, 1)
    m.minimize(s + t)
    with pytest.raises(rs.ModelError) as e:
        m.solve()
    assert e.value.code == "RS4013" and "mip_engine=cpsat" in str(e.value.diagnostics)
    sol = m.solve(options={"mip_engine": "cpsat"})
    assert sol.status == "Optimal" and sol.cp["verified"] is True
    assert sol.objective == pytest.approx(3.0)  # one starts at 0, the other when it is done


def test_cpsat_without_the_library_is_a_typed_error():
    import subprocess
    import sys

    code = (
        "import rusolve as rs\n"
        "m = rs.Model()\n"
        "s = m.add_var(ub=9, integer=True, name='s')\n"
        "m.add_no_overlap([m.add_interval(s, 2)])\n"
        "try:\n"
        "    m.solve(options={'mip_engine': 'cpsat'})\n"
        "except rs.SolveError as e:\n"
        "    print(e.code, str(e))\n"
    )
    out = subprocess.run(
        [sys.executable, "-I", "-c", code],
        env={"RUSOLVE_ORTOOLS_LIB": "/nonexistent/libortools.so.9", "PATH": "/usr/bin:/bin"},
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    assert out.startswith("RS0012")
    assert "/nonexistent/libortools.so.9" in out
    assert "install.sh --with-cpsat" in out and "rusolve[cpsat]" in out


def test_unknown_engine_is_refused():
    m, _ = job_shop()
    with pytest.raises(rs.SolveError) as e:
        m.solve(options={"mip_engine": "nosuch"})
    assert e.value.code == "RS0004"
