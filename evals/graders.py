"""Deterministic graders for the agent-task evals.

Every grader takes the agent's parsed final answer (the JSON after the last `ANSWER:` in its
final message) and returns `(ok, detail)`. Ground truth is recomputed by rusolve itself
(`lib/rusolve.py`): the CLI for solves, `check`, `feasrelax`, and the MCP server for edits. No
expected number is stored in a task file, so the graders stay correct when the solver improves.

`python3 evals/graders.py --selfcheck` proves the graders: a reference answer built from
rusolve's own output must pass and a plausible wrong answer must fail, for every task.
"""
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import rusolve as rs  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def task_dir(task_id):
    return os.path.join(HERE, "tasks", task_id)


def _num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def parse_answer(text):
    """The JSON object after the last `ANSWER:` marker in the agent's final message, or None."""
    if not text or "ANSWER:" not in text:
        return None
    tail = text.rsplit("ANSWER:", 1)[1]
    i = tail.find("{")
    if i < 0:
        return None
    try:
        obj, _ = json.JSONDecoder().raw_decode(tail[i:])
    except ValueError:
        return None
    return obj if isinstance(obj, dict) else None


# ------------------------------------------------------------------ solve_report

def grade_solve_report(ans):
    model = os.path.join(task_dir("solve_report"), "model.mps")
    truth = rs.mip(model)
    for k in ("status", "objective", "selected"):
        if k not in ans:
            return False, f"answer lacks '{k}'"
    if ans["status"] != truth["status"]:
        return False, f"status {ans['status']!r}, rusolve says {truth['status']!r}"
    t = truth["objective"]
    if not _num(ans["objective"]) or abs(ans["objective"] - t) > 1e-6 * (1 + abs(t)):
        return False, f"objective {ans['objective']!r}, rusolve says {t}"
    if not isinstance(ans["selected"], list):
        return False, "'selected' is not a list"
    # The selection must be a feasible solution with that objective (ties make a name-by-name
    # comparison wrong, so it is checked by rusolve's own `check`).
    names = [c["name"] for c in truth["solution"]["columns"]]
    sel = set(ans["selected"])
    unknown = sel - set(names)
    if unknown:
        return False, f"unknown columns in 'selected': {sorted(unknown)}"
    chk = _check(model, {n: 1.0 if n in sel else 0.0 for n in names})
    if chk["verdict"] != "feasible":
        return False, f"selected columns are not feasible: {[f['name'] for f in chk['findings']]}"
    if abs(chk["objective"]["value"] - t) > 1e-6 * (1 + abs(t)):
        return False, f"selected columns have objective {chk['objective']['value']}, not {t}"
    return True, f"status and objective {t} match; selection is feasible"


def _check(model, values):
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump({"columns": [{"name": n, "value": v} for n, v in values.items()]}, f)
        path = f.name
    try:
        return rs.cli_json(["check", model, path, "--format", "json"])
    finally:
        os.unlink(path)


# ------------------------------------------------------------------ debug_infeasible

def grade_debug_infeasible(ans):
    model = os.path.join(task_dir("debug_infeasible"), "model.mps")
    edits = ans.get("edits")
    if not isinstance(edits, list) or not edits:
        return False, "answer needs a non-empty 'edits' list"
    truth = rs.cli_json(["feasrelax", model, "--mode", "min-count", "--format", "json"])
    budget = truth["summary"]["relaxed_count"]
    if len(edits) > budget:
        return False, f"{len(edits)} edits; rusolve's minimum-count repair needs {budget}"
    bad = [e.get("kind") for e in edits if not isinstance(e, dict) or e.get("kind") not in ("col_bounds", "row_bounds")]
    if bad:
        return False, f"only col_bounds / row_bounds edits are allowed, got {bad}"
    with rs.Mcp() as m:
        mid = m.call("load_model", path=model)["model_id"]
        try:
            m.call("modify", model_id=mid, edits=edits)
        except RuntimeError as e:
            return False, f"edits rejected by rusolve: {e}"
        status = m.call("solve", model_id=mid)["status"]
    if status != "Optimal":
        return False, f"model is still {status} after the edits"
    return True, f"{len(edits)} edit(s) make the model Optimal (budget {budget})"


# ------------------------------------------------------------------ verify_formulation

def grade_verify_formulation(ans):
    d = task_dir("verify_formulation")
    truth = rs.cli_json(["check", os.path.join(d, "model.mps"), os.path.join(d, "candidate.json"), "--format", "json"])
    want = sorted({f["name"] for f in truth["findings"]})
    got = ans.get("violations")
    if not isinstance(got, list) or "feasible" not in ans:
        return False, "answer needs 'feasible' and 'violations'"
    if ans["feasible"] != (truth["verdict"] == "feasible"):
        return False, f"feasible={ans['feasible']}, rusolve's verdict is {truth['verdict']}"
    if sorted(set(got)) != want:
        return False, f"violations {sorted(set(got))}, rusolve's check lists {want}"
    return True, f"verdict and violation set {want} match rusolve check"


# ------------------------------------------------------------------ tune_family

TUNE_INSTANCES = ["inst1.mps", "inst2.mps", "inst3.mps"]
TUNE_NODE_FACTOR = 0.8
TUNE_OBJ_TOL = 0.01


def grade_tune_family(ans):
    opts = ans.get("options")
    if not isinstance(opts, dict) or not opts:
        return False, "answer needs a non-empty 'options' object"
    d = task_dir("tune_family")
    base_nodes, new_nodes = 0, 0
    for name in TUNE_INSTANCES:
        path = os.path.join(d, name)
        base = rs.mip(path)
        try:
            run = rs.mip(path, opts)
        except RuntimeError as e:
            return False, f"rusolve rejected the options on {name}: {e}"
        if run["status"] != "Optimal":
            return False, f"{name}: status {run['status']} under the options"
        if abs(run["objective"] - base["objective"]) > TUNE_OBJ_TOL * (1 + abs(base["objective"])):
            return False, f"{name}: objective {run['objective']} is not within 1% of the optimum {base['objective']}"
        base_nodes += base["stats"]["nodes"]
        new_nodes += run["stats"]["nodes"]
    if new_nodes > TUNE_NODE_FACTOR * base_nodes:
        return False, f"total nodes {new_nodes} vs default {base_nodes}: needs at most {TUNE_NODE_FACTOR:.0%} of default"
    return True, f"total nodes {new_nodes} vs default {base_nodes}, all within 1% of optimal"


# ------------------------------------------------------------------ add_cuts

MAX_CUTS = 10
GAP_CLOSED_MIN = 0.5


def _lp_relaxation_file(model):
    """The model with its integrality dropped (it is ours: integers are MARKER blocks and BV bounds)."""
    text = open(model).read().replace(" BV ", " UP ")
    text = "\n".join(l for l in text.splitlines() if "MARKER" not in l) + "\n"
    f = tempfile.NamedTemporaryFile("w", suffix=".mps", delete=False)
    f.write(text)
    f.close()
    return f.name


def grade_add_cuts(ans):
    model = os.path.join(task_dir("add_cuts"), "model.mps")
    cuts = ans.get("cuts")
    if not isinstance(cuts, list) or not cuts or len(cuts) > MAX_CUTS:
        return False, f"answer needs 1..{MAX_CUTS} 'cuts'"
    rows = []
    for i, c in enumerate(cuts):
        if not isinstance(c, dict) or not isinstance(c.get("coefs"), dict) or not c["coefs"]:
            return False, f"cut {i} needs a 'coefs' object {{column: value}}"
        if ("lower" in c and c["lower"] is not None and not _num(c["lower"])) or \
           ("upper" in c and c["upper"] is not None and not _num(c["upper"])):
            return False, f"cut {i}: lower/upper must be numbers or null"
        if c.get("lower") is None and c.get("upper") is None:
            return False, f"cut {i} has neither lower nor upper"
        rows.append({"kind": "add_row", "name": f"agentcut{i}", "lower": c.get("lower"), "upper": c.get("upper"),
                     "coefs": c["coefs"]})
    mip_obj = rs.mip(model)["objective"]
    names = [c["name"] for c in rs.mip(model)["solution"]["columns"]]
    for i, c in enumerate(cuts):
        if set(c["coefs"]) - set(names):
            return False, f"cut {i}: unknown columns {sorted(set(c['coefs']) - set(names))}"
    # Validity: no integer-feasible point of the original model may violate a cut. Optimise the
    # cut's left-hand side over the original MIP (rusolve itself) and compare with its bounds.
    with rs.Mcp() as m:
        mid = m.call("load_model", path=model)["model_id"]
        for i, c in enumerate(cuts):
            costs = lambda s: [{"kind": "cost", "name": n, "cost": s * float(c["coefs"].get(n, 0.0))} for n in names]
            if c.get("upper") is not None:
                m.call("modify", model_id=mid, edits=costs(-1.0))  # minimise -a.x = maximise a.x
                best = -m.call("solve", model_id=mid)["objective"]
                if best > c["upper"] + 1e-6:
                    return False, f"cut {i} is invalid: a feasible point reaches {best} > upper {c['upper']}"
            if c.get("lower") is not None:
                m.call("modify", model_id=mid, edits=costs(1.0))
                best = m.call("solve", model_id=mid)["objective"]
                if best < c["lower"] - 1e-6:
                    return False, f"cut {i} is invalid: a feasible point reaches {best} < lower {c['lower']}"
    lp_file = _lp_relaxation_file(model)
    try:
        with rs.Mcp() as m:
            mid = m.call("load_model", path=lp_file)["model_id"]
            lp0 = m.call("solve", model_id=mid)["objective"]
            m.call("modify", model_id=mid, edits=rows)
            lp1 = m.call("solve", model_id=mid)["objective"]
    finally:
        os.unlink(lp_file)
    gap = mip_obj - lp0
    closed = (lp1 - lp0) / gap if gap > 1e-9 else 1.0
    if closed < GAP_CLOSED_MIN:
        return False, f"cuts close {closed:.0%} of the root gap (LP {lp0:.4f} -> {lp1:.4f}, MIP {mip_obj:.4f}); need {GAP_CLOSED_MIN:.0%}"
    return True, f"{len(cuts)} valid cut(s) close {closed:.0%} of the root gap"


# ------------------------------------------------------------------ extend_wasm

# Both solves use the HiGHS child rule (mip_child_selection=-1): the default hybrid rule dives so
# much that the open-node queue is consulted only about a dozen times, which leaves a node-selection
# rule almost nothing to decide. Under the HiGHS rule the reference (selfcheck/reference_depth_first.wat)
# reaches 534 nodes against 912 for the best-bound pop, so 0.8 leaves room; a component with no
# opinion leaves the tree at exactly 912 and fails. Measured numbers are in the task's section of
# README.md.
EXTEND_NODE_FACTOR = 0.8
EXTEND_OBJ_TOL = 1e-6
EXTEND_OPTIONS = {"mip_child_selection": -1}


def grade_extend_wasm(ans):
    model = os.path.join(task_dir("extend_wasm"), "model.mps")
    path = ans.get("component_path")
    if not isinstance(path, str) or not path.strip():
        return False, "answer needs 'component_path', the path of the component it wrote"
    if not os.path.isfile(path):
        return False, f"no component file at {path!r}"
    base = rs.mip(model, EXTEND_OPTIONS)
    try:
        run = rs.mip(model, EXTEND_OPTIONS, wasm_hook=path)  # a component rusolve cannot load exits 3 -> failure
    except RuntimeError as e:
        return False, f"rusolve would not run the component: {e}"
    if run["status"] != "Optimal":
        return False, f"with the component the model is {run['status']}, not Optimal"
    t = base["objective"]
    if abs(run["objective"] - t) > EXTEND_OBJ_TOL * (1 + abs(t)):
        return False, (f"the component changed the answer: objective {run['objective']} vs "
                       f"{t} without it (a hook may only reorder the search)")
    w = run["stats"].get("wasm")
    if not isinstance(w, dict):
        return False, "the result has no stats.wasm: the component was not used"
    if not w.get("calls"):
        return False, f"the component was never called: {w}"
    for k in ("traps", "out_of_fuel", "invalid"):
        if w.get(k):
            return False, f"the component {k} {w[k]} times: {w}"
    if not w.get("ok"):
        return False, f"the component never made a decision (declined {w.get('declined')}): {w}"
    n0, n1 = base["stats"]["nodes"], run["stats"]["nodes"]
    if n1 > EXTEND_NODE_FACTOR * n0:
        return False, (f"{n1} nodes with the component vs {n0} without: needs at most "
                       f"{EXTEND_NODE_FACTOR:.0%} of the solve without it")
    return True, (f"{n1} nodes vs {n0} without it, {w['ok']}/{w['calls']} calls answered, "
                  f"objective {t} unchanged")


# ------------------------------------------------------------------ class_* (playbook skills)
# One task per problem-class playbook (docs/skills/<class>). Ground truth is recomputed by rusolve
# (`mip` envelope, `lint`); the playbook name is the one fixed fact: it is the skill that class maps to.

PLAYBOOK = {
    "class_flow": "rusolve-assignment-transport-flow", "class_cover": "rusolve-set-cover-packing",
    "class_scheduling": "rusolve-scheduling", "class_lotsizing": "rusolve-lot-sizing",
    "class_facility": "rusolve-facility-location", "class_routing": "rusolve-routing",
    "class_pwl_semi": "rusolve-piecewise-semi", "class_qp": "rusolve-convex-qp",
    "class_conic": "rusolve-conic", "class_miqp": "rusolve-convex-qp", "class_misocp": "rusolve-conic",
}
# The two mixed-integer classes are refused by default; their experimental engine is a lever, and the
# truth is what rusolve answers with the lever on (the refusal message names it).
LEVER = {"class_miqp": {"miqp_oa": True}, "class_misocp": {"misocp_oa": True}}
TSP_PTS = [(0, 0), (4, 1), (7, 5), (3, 8), (-2, 5), (-3, 1)]   # gen_models.class_routing


def _class_truth(tid):
    return rs.mip(os.path.join(task_dir(tid), "model.lp"), options=LEVER.get(tid))


def _class_common(tid, ans, truth, tol=1e-6):
    """status, objective and playbook agree with rusolve; returns an error string or None."""
    for k in ("status", "objective", "playbook"):
        if k not in ans:
            return f"answer lacks '{k}'"
    if ans["status"] != truth["status"]:
        return f"status {ans['status']!r}, rusolve says {truth['status']!r}"
    t = truth["objective"]
    if not _num(ans["objective"]) or abs(ans["objective"] - t) > tol * (1 + abs(t)):
        return f"objective {ans['objective']!r}, rusolve says {t}"
    if ans["playbook"] != PLAYBOOK[tid]:
        return f"playbook {ans['playbook']!r}, this problem class is covered by {PLAYBOOK[tid]!r}"
    return None


def _class_field(ans, truth, *keys):
    for k in keys:
        if ans.get(k) != truth.get(k):
            return f"{k} {ans.get(k)!r}, rusolve's envelope says {truth.get(k)!r}"
    return None


def _grade_class_plain(tid, ans):
    truth = _class_truth(tid)
    err = _class_common(tid, ans, truth) or (_class_field(ans, truth, "class") if "class" in ans or tid in ("class_flow", "class_cover") else None)
    return (False, err) if err else (True, f"{truth['status']} {truth['objective']}, class {truth['class']}, playbook {PLAYBOOK[tid]}")


def grade_class_flow(ans):
    return _grade_class_plain("class_flow", ans)


def grade_class_cover(ans):
    return _grade_class_plain("class_cover", ans)


def _lint_bigm(tid):
    lint = rs.cli_json(["lint", os.path.join(task_dir(tid), "model.lp"), "--format", "json"])
    return sum(f.get("count", 0) for f in lint["findings"] if f["rule"] == "big-m")


def _grade_class_bigm(tid, ans):
    truth = _class_truth(tid)
    err = _class_common(tid, ans, truth)
    if err:
        return False, err
    want = _lint_bigm(tid)
    if ans.get("big_m_rows") != want:
        return False, f"big_m_rows {ans.get('big_m_rows')!r}, rusolve's lint flags {want} rows"
    return True, f"{truth['objective']} and {want} big-M rows flagged, as rusolve lint reports"


def grade_class_scheduling(ans):
    return _grade_class_bigm("class_scheduling", ans)


def grade_class_lotsizing(ans):
    return _grade_class_bigm("class_lotsizing", ans)


def grade_class_facility(ans):
    truth = _class_truth("class_facility")
    err = _class_common("class_facility", ans, truth)
    if err:
        return False, err
    want = sorted(c["name"] for c in truth["solution"]["columns"] if c["name"].startswith("open_") and c["value"] > 0.5)
    got = ans.get("open")
    if not isinstance(got, list) or sorted(got) != want:
        return False, f"open {got!r}, rusolve opens {want}"
    return True, f"objective {truth['objective']}, open sites {want}"


def grade_class_routing(ans):
    truth = _class_truth("class_routing")
    err = _class_common("class_routing", ans, truth)
    if err:
        return False, err
    t = ans.get("tour")
    n = len(TSP_PTS)
    if not isinstance(t, list) or sorted(t) != list(range(n)) or t[0] != 0:
        return False, f"tour {t!r} is not one cycle through nodes 0..{n - 1} starting at 0"
    d = lambda a, b: round(((TSP_PTS[a][0] - TSP_PTS[b][0]) ** 2 + (TSP_PTS[a][1] - TSP_PTS[b][1]) ** 2) ** 0.5 * 10)
    length = sum(d(a, b) for a, b in zip(t, t[1:] + [0]))
    if abs(length - truth["objective"]) > 1e-6:
        return False, f"the tour has length {length}, the optimum is {truth['objective']}"
    return True, f"tour {t} has the optimal length {length}"


def grade_class_pwl_semi(ans):
    truth = _class_truth("class_pwl_semi")
    err = _class_common("class_pwl_semi", ans, truth) or _class_field(ans, truth, "class", "engine")
    return (False, err) if err else (True, f"{truth['status']} {truth['objective']} via {truth['class']}/{truth['engine']}")


def grade_class_qp(ans):
    truth = _class_truth("class_qp")
    err = _class_common("class_qp", ans, truth, tol=1e-4) or _class_field(ans, truth, "class")
    if err:
        return False, err
    want = {c["name"]: c["value"] for c in truth["solution"]["columns"]}
    w = ans.get("weights")
    if not isinstance(w, list) or len(w) != 4 or not all(_num(v) for v in w):
        return False, "'weights' must be a list of four numbers"
    for i, v in enumerate(w):
        if abs(v - want[f"w{i}"]) > 2e-3:
            return False, f"w{i} = {v}, rusolve's solution has {want[f'w{i}']:.4f}"
    return True, f"{truth['class']} objective {truth['objective']:.6g} and weights match"


def _values(ans, truth, names, tol):
    """`values` (a name -> number map) matches rusolve's solution on `names`; an error string or None."""
    want = {c["name"]: c["value"] for c in truth["solution"]["columns"]}
    v = ans.get("values")
    if not isinstance(v, dict) or not all(_num(v.get(n)) for n in names):
        return f"'values' must map {', '.join(names)} to numbers"
    for n in names:
        if abs(v[n] - want[n]) > tol:
            return f"{n} = {v[n]}, rusolve's solution has {want[n]:.4f}"
    return None


def _grade_class_values(tid, ans):
    truth = _class_truth(tid)
    err = (_class_common(tid, ans, truth, tol=1e-4) or _class_field(ans, truth, "class", "engine")
           or _values(ans, truth, "abc", 2e-3))
    return (False, err) if err else (True, f"{truth['class']}/{truth['engine']} objective {truth['objective']:.6g} and values match")


def grade_class_conic(ans):
    return _grade_class_values("class_conic", ans)


def grade_class_misocp(ans):
    return _grade_class_values("class_misocp", ans)


def grade_class_miqp(ans):
    truth = _class_truth("class_miqp")
    err = _class_common("class_miqp", ans, truth, tol=1e-4) or _class_field(ans, truth, "class", "engine")
    if err:
        return False, err
    cols = {c["name"]: c["value"] for c in truth["solution"]["columns"]}
    w = ans.get("weights")
    if not isinstance(w, list) or len(w) != 4 or not all(_num(v) for v in w):
        return False, "'weights' must be a list of four numbers"
    for i, v in enumerate(w):
        if abs(v - cols[f"w{i}"]) > 2e-3:
            return False, f"w{i} = {v}, rusolve's solution has {cols[f'w{i}']:.4f}"
    held = sorted(n for n, v in cols.items() if n.startswith("y") and v > 0.5)
    if sorted(ans.get("held") or []) != held:
        return False, f"held {ans.get('held')!r}, rusolve's solution holds {held}"
    return True, f"{truth['class']}/{truth['engine']} objective {truth['objective']:.6g}, held {held}"


ROUTER_MODELS = {"m1": "class_flow", "m2": "class_cover", "m3": "class_qp", "m4": "class_pwl_semi"}


def grade_class_router(ans):
    models = ans.get("models")
    if not isinstance(models, dict):
        return False, "answer needs a 'models' object with m1..m4"
    for m, src in ROUTER_MODELS.items():
        got = models.get(m)
        if not isinstance(got, dict):
            return False, f"no entry for {m}"
        truth = rs.mip(os.path.join(task_dir("class_router"), f"{m}.lp"))
        if got.get("class") != truth["class"]:
            return False, f"{m}: class {got.get('class')!r}, rusolve says {truth['class']!r}"
        if got.get("playbook") != PLAYBOOK[src]:
            return False, f"{m}: playbook {got.get('playbook')!r}, expected {PLAYBOOK[src]!r}"
    return True, "all four models routed to the right class and playbook"


# ------------------------------------------------------------------ word problems
# The agent models a word problem in Python (`rusolve.Model`). Ground truth: the reference model
# in tasks/<id>/selfcheck/reference.mps (the same problem, written by hand), solved by rusolve at
# grading time. The agent's reported plan is then judged against that model with rusolve's own
# `check`, so a plan that breaks a stated rule fails even when its objective looks right.

WORD_PROBLEMS = ("diet_plan", "project_selection", "fixed_charge_supply", "shift_cover", "mcp_build_model",
                 "supplier_sourcing")


def _reference(tid):
    """The hand-written reference model: MPS, or LP when the problem has general constraints."""
    for ext in ("mps", "lp"):
        p = os.path.join(task_dir(tid), "selfcheck", f"reference.{ext}")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"tasks/{tid}/selfcheck/reference.(mps|lp)")


def _grade_word_problem(tid, ans):
    model = _reference(tid)
    truth = rs.mip(model)
    for k in ("status", "objective", "plan"):
        if k not in ans:
            return False, f"answer lacks '{k}'"
    if ans["status"] != truth["status"]:
        return False, f"status {ans['status']!r}, rusolve says {truth['status']!r}"
    t = truth["objective"]
    if not _num(ans["objective"]) or abs(ans["objective"] - t) > 1e-4 * (1 + abs(t)):
        return False, f"objective {ans['objective']!r}, rusolve's reference model says {t}"
    plan = ans["plan"]
    if not isinstance(plan, dict) or not all(_num(v) for v in plan.values()):
        return False, "'plan' must map variable names to numbers"
    names = [c["name"] for c in truth["solution"]["columns"]]
    unknown = sorted(set(plan) - set(names))
    if unknown:
        return False, f"unknown variables in 'plan': {unknown}"
    chk = _check(model, {n: float(plan.get(n, 0.0)) for n in names})
    if chk["verdict"] != "feasible":
        return False, f"plan breaks the problem's rules: {sorted({f['name'] for f in chk['findings']})}"
    if abs(chk["objective"]["value"] - t) > 1e-4 * (1 + abs(t)):
        return False, f"plan has objective {chk['objective']['value']}, the optimum is {t}"
    return True, f"status and objective {t:.6g} match; the plan is feasible and optimal"


def grade_diet_plan(ans):
    return _grade_word_problem("diet_plan", ans)


def grade_project_selection(ans):
    return _grade_word_problem("project_selection", ans)


def grade_fixed_charge_supply(ans):
    return _grade_word_problem("fixed_charge_supply", ans)


def grade_shift_cover(ans):
    return _grade_word_problem("shift_cover", ans)


def grade_mcp_build_model(ans):
    return _grade_word_problem("mcp_build_model", ans)


def grade_supplier_sourcing(ans):
    return _grade_word_problem("supplier_sourcing", ans)


# ------------------------------------------------------------------ class_cp_scheduling
# Native scheduling objects have no file syntax, so the reference is built in Python by
# `cp_reference.py` (run with RSV_PY, an interpreter that has the `rusolve` package), solved at
# grading time, and the reported plan is judged with `Model.check`, which evaluates the intervals
# and no-overlap objects themselves.

def _cp_ref(*args):
    import subprocess
    py = os.environ.get("RSV_PY", "python3")
    r = subprocess.run([py, "-I", os.path.join(HERE, "cp_reference.py"), *args], capture_output=True, text=True,
                       timeout=300)
    if r.returncode != 0:
        raise RuntimeError(f"cp_reference.py {args[0]} exited {r.returncode}: {r.stderr.strip()[-400:]}")
    return json.loads(r.stdout)


def grade_class_cp_scheduling(ans):
    truth = _cp_ref("solve")
    for k in ("status", "objective", "plan"):
        if k not in ans:
            return False, f"answer lacks '{k}'"
    if ans["status"] != truth["status"]:
        return False, f"status {ans['status']!r}, rusolve says {truth['status']!r}"
    t = truth["objective"]
    if not _num(ans["objective"]) or abs(ans["objective"] - t) > 1e-6 * (1 + abs(t)):
        return False, f"objective {ans['objective']!r}, rusolve's reference model says {t}"
    plan = ans["plan"]
    if not isinstance(plan, dict) or not all(_num(v) for v in plan.values()):
        return False, "'plan' must map variable names to numbers"
    reported = {n for n in truth["plan"] if n.startswith("s") or n == "cmax"}
    unknown = sorted(set(plan) - reported)
    if unknown:
        return False, f"unknown variables in 'plan': {unknown}"
    chk = _cp_ref("check", json.dumps(plan))
    if chk["verdict"] != "feasible":
        return False, f"plan breaks the problem's rules: {chk['findings']}"
    if abs(chk["objective"] - t) > 1e-6 * (1 + abs(t)):
        return False, f"plan has objective {chk['objective']}, the optimum is {t}"
    return True, f"status and objective {t:.6g} match; the schedule is feasible on the intervals and no-overlap objects"


GRADERS = {
    "solve_report": grade_solve_report,
    "debug_infeasible": grade_debug_infeasible,
    "verify_formulation": grade_verify_formulation,
    "tune_family": grade_tune_family,
    "add_cuts": grade_add_cuts,
    "extend_wasm": grade_extend_wasm,
    "class_router": grade_class_router,
    "class_flow": grade_class_flow,
    "class_cover": grade_class_cover,
    "class_scheduling": grade_class_scheduling,
    "class_lotsizing": grade_class_lotsizing,
    "class_facility": grade_class_facility,
    "class_routing": grade_class_routing,
    "class_pwl_semi": grade_class_pwl_semi,
    "class_qp": grade_class_qp,
    "class_conic": grade_class_conic,
    "class_miqp": grade_class_miqp,
    "class_misocp": grade_class_misocp,
    "diet_plan": grade_diet_plan,
    "project_selection": grade_project_selection,
    "fixed_charge_supply": grade_fixed_charge_supply,
    "shift_cover": grade_shift_cover,
    "mcp_build_model": grade_mcp_build_model,
    "supplier_sourcing": grade_supplier_sourcing,
    "class_cp_scheduling": grade_class_cp_scheduling,
}


def grade(task_id, answer_text):
    ans = parse_answer(answer_text)
    if ans is None:
        return False, "no parsable 'ANSWER: {json}' in the final message"
    try:
        return GRADERS[task_id](ans)
    except Exception as e:  # a grader bug must never read as a pass
        return False, f"grader error: {type(e).__name__}: {e}"


# ------------------------------------------------------------------ selfcheck

def _class_cases():
    """Reference (built from rusolve's envelope and lint) and a wrong answer for each class_* task."""
    out = {}
    for tid in PLAYBOOK:
        t = _class_truth(tid)
        base = {"status": t["status"], "objective": t["objective"], "playbook": PLAYBOOK[tid]}
        wrong = dict(base, objective=t["objective"] + 1)
        good = dict(base)
        if tid in ("class_flow", "class_cover", "class_pwl_semi", "class_qp", "class_conic", "class_miqp", "class_misocp"):
            good["class"] = t["class"]
        if tid in ("class_pwl_semi", "class_conic", "class_miqp", "class_misocp"):
            good["engine"] = t["engine"]
        if tid in ("class_conic", "class_misocp"):
            cols = {c["name"]: c["value"] for c in t["solution"]["columns"]}
            good["values"] = {n: round(cols[n], 4) for n in "abc"}
            wrong = dict(good, values=dict(good["values"], a=good["values"]["a"] + 0.5))
        if tid == "class_miqp":
            cols = {c["name"]: c["value"] for c in t["solution"]["columns"]}
            good["weights"] = [round(cols[f"w{i}"], 4) for i in range(4)]
            good["held"] = sorted(n for n, v in cols.items() if n.startswith("y") and v > 0.5)
            wrong = dict(good, held=["y0", "y1", "y2", "y3"])
        if tid in ("class_scheduling", "class_lotsizing"):
            good["big_m_rows"] = _lint_bigm(tid)
            wrong = dict(good, big_m_rows=0)
        if tid == "class_facility":
            good["open"] = sorted(c["name"] for c in t["solution"]["columns"]
                                  if c["name"].startswith("open_") and c["value"] > 0.5)
            wrong = dict(good, playbook="rusolve-routing")
        if tid == "class_routing":
            nxt = {}
            for c in t["solution"]["columns"]:
                if c["name"].startswith("x_") and c["value"] > 0.5:
                    _, a, b = c["name"].split("_")
                    nxt[int(a)] = int(b)
            tour, k = [0], nxt[0]
            while k != 0:
                tour.append(k)
                k = nxt[k]
            good["tour"] = tour
            wrong = dict(good, tour=tour[:-1])   # misses a node: not a tour
        if tid == "class_qp":
            cols = {c["name"]: c["value"] for c in t["solution"]["columns"]}
            good["weights"] = [round(cols[f"w{i}"], 4) for i in range(4)]
            wrong = dict(good, weights=[0.25, 0.25, 0.25, 0.25])
        out[tid] = (good, wrong)
    routed = {}
    for m, src in ROUTER_MODELS.items():
        routed[m] = {"class": rs.mip(os.path.join(task_dir("class_router"), f"{m}.lp"))["class"], "playbook": PLAYBOOK[src]}
    bad = {k: dict(v) for k, v in routed.items()}
    bad["m4"]["playbook"] = "rusolve-convex-qp"
    out["class_router"] = ({"models": routed}, {"models": bad})
    return out


def _rule_breakers(tid, sol):
    """(label, plan, row) wrong plans whose objective equals the optimum but that break a structural
    row; the grader must reject each and name that row. Found by hand, checked at run time."""
    if tid == "diet_plan":   # 0.6 beans for 0.5 bread: same cost and calories, beans above 40% of calories
        return [("beans_share broken", {**sol, "beans": sol["beans"] + 0.5, "bread": sol["bread"] - 0.6}, "beans_share")]
    if tid == "project_selection":
        return [("E_needs_A broken", {"pB": 1, "pC": 1, "pE": 1, "pF": 2 / 9}, "E_needs_A"),   # E without A
                ("C_or_D broken", {"pA": 1, "pB": 1, "pC": 1, "pD": 1}, "C_or_D")]
    if tid == "fixed_charge_supply":   # 5 units for store 2 from closed warehouse 2, repaid on store 4
        return [("cap_2 broken", {**sol, "ship_3_2": 20, "ship_2_2": 5, "ship_3_4": 5, "ship_1_4": 30}, "cap_2")]
    if tid == "shift_cover":   # one Monday starter moved to Tuesday: same cost, Monday uncovered
        return [("cover_mon broken", {**sol, "start_mon": sol["start_mon"] - 1, "start_tue": sol["start_tue"] + 1}, "cover_mon")]
    if tid == "supplier_sourcing":   # B delivers with its switch off; A and C both deliver; B below its minimum lot
        return [("b_off broken", {**sol, "useB": 0}, "b_off"),
                ("a_or_c broken", {"qA": 10, "qB": 70, "qC": 20, "useB": 1, "largest": 70}, "a_or_c"),
                ("b_min broken", {"qB": 10, "qC": 50, "qA": 0, "useB": 1, "largest": 50}, "b_min")]
    if tid == "mcp_build_model":   # 17 rolls on 5 loaves: the same profit, rolls above 1.5 x loaves
        return [("rolls_cap broken", {"loaves": 5, "rolls": 17, "cakes": 1}, "rolls_cap")]
    return []


def _selfcheck():
    """Reference answers (built from rusolve) must pass; plausible wrong ones must fail."""
    d = task_dir
    t = rs.mip(os.path.join(d("solve_report"), "model.mps"))
    sel = [c["name"] for c in t["solution"]["columns"] if c["value"] > 0.5]
    fe = rs.cli_json(["feasrelax", os.path.join(d("debug_infeasible"), "model.mps"), "--mode", "min-count",
                      "--format", "json"])["edits"][0]
    ref_edit = {"kind": "row_bounds", "name": fe["name"], "lower": fe["to"]}
    ck = rs.cli_json(["check", os.path.join(d("verify_formulation"), "model.mps"),
                      os.path.join(d("verify_formulation"), "candidate.json"), "--format", "json"])
    names = [c["name"] for c in rs.mip(os.path.join(d("add_cuts"), "model.mps"))["solution"]["columns"]]
    cases = {
        "solve_report": (
            {"status": t["status"], "objective": t["objective"], "selected": sel},
            {"status": "Optimal", "objective": t["objective"] + 1, "selected": sel}),
        "debug_infeasible": (
            {"edits": [ref_edit]},
            {"edits": [{"kind": "row_bounds", "name": "cap1", "upper": 50}]}),
        "verify_formulation": (
            {"feasible": False, "violations": sorted({f["name"] for f in ck["findings"]})},
            {"feasible": False, "violations": ["cap"]}),
        "tune_family": (
            {"options": {"mip_rel_gap": 0.01}},
            {"options": {"mip_heuristic_effort": 0.0}}),
        "add_cuts": (
            {"cuts": [{"name": "k", "upper": 3, "coefs": {n: 1 for n in names}}]},
            {"cuts": [{"name": "k", "upper": 1, "coefs": {n: 1 for n in names}}]}),  # invalid: cuts off optimum
        # The agent's code runs inside the solver (Wasm hook). Fixtures under tasks/extend_wasm/selfcheck.
        "extend_wasm": (
            {"component_path": os.path.join(d("extend_wasm"), "selfcheck", "reference_depth_first.wat")},
            [("no opinion", {"component_path": os.path.join(d("extend_wasm"), "selfcheck", "no_opinion.wat")}),
             ("traps", {"component_path": os.path.join(d("extend_wasm"), "selfcheck", "traps.wat")}),
             ("not loadable", {"component_path": os.path.join(d("extend_wasm"), "selfcheck",
                                                              "imports_host_function.wat")}),
             ("no component", {"component_path": "does-not-exist.wat"})]),
    }
    cases.update(_class_cases())
    for tid in WORD_PROBLEMS:
        ref = rs.mip(_reference(tid))
        sol = {c["name"]: c["value"] for c in ref["solution"]["columns"] if abs(c["value"]) > 1e-9}
        first = next(iter(sol))
        cases[tid] = (
            {"status": ref["status"], "objective": ref["objective"], "plan": sol},
            [("objective off", {"status": ref["status"], "objective": ref["objective"] * 1.1 + 1, "plan": sol}),
             ("plan breaks a rule", {"status": ref["status"], "objective": ref["objective"],
                                     "plan": {**sol, first: sol[first] + 1000}}),
             ("no plan", {"status": ref["status"], "objective": ref["objective"]}),
             *[(label, {"status": ref["status"], "objective": ref["objective"], "plan": plan}, row)
               for label, plan, row in _rule_breakers(tid, sol)]])
    ref = _cp_ref("solve")
    sol = {n: v for n, v in ref["plan"].items() if n.startswith("s") or n == "cmax"}
    head = {"status": ref["status"], "objective": ref["objective"]}

    def pulled_to_first(k):
        """The operation of machine k that starts last, moved to the earliest start on that machine:
        the same makespan (it only moves earlier), but it now overlaps another operation there."""
        ops = sorted((n for n in sol if n.startswith("s") and n.endswith(str(k))), key=sol.get)
        return {**sol, ops[-1]: sol[ops[0]]}

    cases["class_cp_scheduling"] = (
        {**head, "plan": sol},
        [("objective off", {**head, "objective": ref["objective"] + 1, "plan": sol}),
         ("no plan", head),
         ("overlap on M1", {**head, "plan": pulled_to_first(1)}, "machine1"),
         ("overlap on M2", {**head, "plan": pulled_to_first(2)}, "machine2")])
    failed = 0
    for tid, (good, bad) in cases.items():
        g_ok, g_detail = GRADERS[tid](good)
        bads = bad if isinstance(bad, list) else [("wrong answer", bad)]
        results = [(item[0], *GRADERS[tid](item[1])) for item in bads]
        # a rule-breaker (third element) must be rejected for that row, not for some other reason
        named = all(len(item) < 3 or item[2] in why for item, (_, _, why) in zip(bads, results))
        verdict = "ok" if g_ok and named and not any(ok for _, ok, _ in results) else "BROKEN"
        failed += verdict != "ok"
        detail = " | ".join(f"{label}: {why}" for label, _, why in results)
        print(f"{tid:20s} {verdict}  reference: {g_detail} | {detail}")
    return failed


if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        sys.exit(1 if _selfcheck() else 0)
    print(__doc__)
