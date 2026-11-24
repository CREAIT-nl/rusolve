"""Reference model and judge of the `class_cp_scheduling` task.

The task states a small job shop; the agent builds it with the native scheduling objects of the
Python API (`add_interval`, `add_no_overlap`). The scheduling objects have no file syntax, so the
reference cannot be a file: this script builds the same problem with the `rusolve` Python package
(run it with an interpreter that has the package, `RSV_PY`), solves it, and judges a reported plan
with `Model.check`, which evaluates the intervals and no-overlap objects themselves.

    python3 cp_reference.py solve                 -> {"status", "objective", "plan"}
    python3 cp_reference.py check '<plan json>'   -> {"verdict", "objective", "findings"}

Not importable by the graders (their `rusolve` is the CLI wrapper in `lib/`); they call it as a
subprocess.
"""
import json
import sys

import rusolve as rs

# job -> (time on machine 1, time on machine 2); machine 1 first, then machine 2
JOBS = {"A": (2, 3), "B": (3, 1), "C": (2, 2)}
HORIZON = sum(a + b for a, b in JOBS.values())


def build():
    """The reference model. Columns: the starts `s<job><machine>` and `cmax` the prompt names, and
    the ends `e<job><machine>` that the reference needs and a plan does not report."""
    m = rs.Model(name="job_shop")
    cmax = m.add_var(ub=HORIZON, integer=True, cost=1, name="cmax")
    ivs = {}
    for job, (d1, d2) in JOBS.items():
        for k, d in ((1, d1), (2, d2)):
            s = m.add_var(ub=HORIZON, integer=True, name=f"s{job}{k}")
            e = m.add_var(ub=HORIZON, integer=True, name=f"e{job}{k}")
            ivs[job, k] = m.add_interval(s, d, end=e, name=f"{job}{k}")
        m.add_constr(ivs[job, 2].start >= ivs[job, 1].end)
        m.add_constr(cmax >= ivs[job, 2].end)
    m.add_no_overlap([ivs[j, 1] for j in JOBS], name="machine1")
    m.add_no_overlap([ivs[j, 2] for j in JOBS], name="machine2")
    return m


def values_of(m, plan):
    """A full assignment from a plan of starts and `cmax`: each end is its start plus the size."""
    size = {f"{j}{k}": JOBS[j][k - 1] for j in JOBS for k in (1, 2)}
    x = {}
    for j in range(m.num_vars):
        name = m.var_name(m.var(j))
        if name.startswith("e"):
            x[name] = float(plan.get("s" + name[1:], 0.0)) + size[name[1:]]
        else:
            x[name] = float(plan.get(name, 0.0))
    return x


def main():
    m = build()
    if sys.argv[1] == "solve":
        sol = m.solve()
        assert sol.cp is not None and sol.cp["verified"]
        plan = {m.var_name(m.var(j)): float(sol.x[j]) for j in range(m.num_vars)}
        print(json.dumps({"status": sol.status, "objective": sol.objective, "plan": plan}))
        return
    plan = json.loads(sys.argv[2])
    r = m.check(values_of(m, plan))
    print(json.dumps({"verdict": r["verdict"], "objective": r["objective"]["value"],
                      "findings": sorted({f["name"] for f in r["findings"]})}))


if __name__ == "__main__":
    main()
