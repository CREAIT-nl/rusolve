#!/usr/bin/env python3
"""MIP correctness gate: `rusolve-cli mip` vs HiGHS (cross-checked by SCIP).

Case set (written to a fresh tools/diff/mip_cases/run-<time>-<pid>/ per run,
gitignored, or --cases-dir DIR):
  (a) seeded random small MIPs (n <= 15, m <= 12): random integer subsets
      (binaries, general ints with finite / infinite bounds, free ints),
      equality + ranged rows, maximize, plus constructed infeasible-by-parity,
      LP-unbounded-but-MIP-infeasible and unbounded-MIP cases;
  (b) the `small` size (generators.MIP_SMALL) of every MIP generator family
      and its variants x a few seeds.

Oracle: highspy (mip_rel_gap=0, tiny abs gap, threads=1, fixed seed). SCIP
(pyscipopt) solves every case too; a HiGHS/SCIP mismatch is reported as
ORACLE-DISAGREE, never as our failure — unless HiGHS with presolve OFF
agrees with SCIP, in which case that 2-of-3 majority becomes the reference
(flagged HIGHS-PRESOLVE-BUG; HiGHS 1.15.1 presolve has given wrong statuses
on rare random cases). Every oracle solution and every
generator's planted witness is re-verified by the independent checker in
mip_check.py, which re-reads the MPS file itself.

Pass = status matches HiGHS (Infeasible/Unbounded both accepted only when
HiGHS reports kUnboundedOrInfeasible) AND, for Optimal:
|obj - hobj| <= 1e-6 max(1, |hobj|), x passes the checker, and the
checker-recomputed objective matches our reported one.

Our binary: $RUSOLVE_BIN (default target/release/rusolve-cli), contract
    rusolve-cli mip <file.mps> --csv [--node-limit N] [--time-limit S]
    -> header `status,objective,bound,nodes,lp_iters,time_s,x`, then one row.

Usage: python tools/diff/mip_diff.py [--oracle-only] [--cases N]
       [--seed S] [--family NAME] [--only-random] [--timeout SEC]
Our binary gets --time-limit 0.8*timeout. A limit hit on an Infeasible /
InfOrUnb case with an integer column of infinite bound is reported as
EXPECTED-NONTERMINATION (plain B&B cannot prove it without gcd reasoning)
and does not fail the gate; every other limit hit is a FAIL.
Exit 0 iff no FAIL and no ORACLE-DISAGREE.
"""
from __future__ import annotations

import argparse
import math
import multiprocessing
import os
import random
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CASES = os.path.join(HERE, "mip_cases")
sys.path.insert(0, HERE)

import highspy  # noqa: E402
import pyscipopt  # noqa: E402
from highspy import HighsModelStatus as HMS  # noqa: E402

from harness import generators as G  # noqa: E402
from mip_check import INF, check_solution, parse_vector, read_mps, row_bounds, write_mps  # noqa: E402

ORACLE_TIME = 120.0
IORU = "InfOrUnb"


# --------------------------------------------------------------------------- cases


def _empty_case(name, n, m):
    return dict(name=name, n=n, m=m, rows=[[] for _ in range(m)], senses=["L"] * m,
                rhs=[0.0] * m, ranges=[0.0] * m, obj=[0.0] * n, lb=[0.0] * n, ub=[INF] * n,
                is_int=[False] * n, maximize=False, bv=False, family="random",
                planted=None, expect=None)


def gen_random_case(rng: random.Random, k: int) -> dict:
    """Random small MIP; ~15% of cases are constructed special statuses."""
    kind = rng.random()
    if kind < 0.06:
        return _parity_infeasible(rng, k, lp_unbounded=False)
    if kind < 0.10:
        return _parity_infeasible(rng, k, lp_unbounded=True)
    if kind < 0.14:
        return _unbounded_mip(rng, k)
    n, m = rng.randint(3, 15), rng.randint(2, 12)
    c = _empty_case(f"rand_{k:03d}", n, m)
    c["maximize"] = rng.random() < 0.3
    c["bv"] = rng.random() < 0.5
    for j in range(n):
        c["obj"][j] = float(rng.choice([-3, -2, -1, -0.5, 0, 1, 1.5, 2, 3, 5]))
        integer = rng.random() < 0.6
        c["is_int"][j] = integer
        r = rng.random()
        if integer and r < 0.4:
            c["lb"][j], c["ub"][j] = 0.0, 1.0
        elif integer and r < 0.7:
            lo = float(rng.choice([0, 0, -2, 1, -5]))
            c["lb"][j], c["ub"][j] = lo, lo + float(rng.choice([1, 3, 5, 10]))
        elif integer and r < 0.85:
            c["lb"][j], c["ub"][j] = float(rng.choice([0, 2, -3])), INF
        elif integer and r < 0.93:
            c["lb"][j], c["ub"][j] = -INF, INF
        elif integer:
            c["lb"][j], c["ub"][j] = -INF, float(rng.choice([0, 4, 8]))
        else:
            c["lb"][j] = float(rng.choice([0, 0, 0, 1, 2.5, -INF]))
            c["ub"][j] = INF if rng.random() < 0.4 else (
                (c["lb"][j] if c["lb"][j] > -INF else 0.0) + float(rng.choice([2, 4.5, 10])))
        if rng.random() < 0.04:  # fixed column
            v = c["lb"][j] if c["lb"][j] > -INF else 1.0
            c["lb"][j] = c["ub"][j] = v
    # Planted point (92% of cases): integer on integer columns, multiples of
    # 0.5 on continuous ones; rhs/ranges are then set around its activity so
    # the case is feasible by construction. The other 8% keep random rhs
    # (mostly infeasible, unstructured).
    planted = rng.random() < 0.92
    x0 = []
    for j in range(n):
        lo = c["lb"][j] if c["lb"][j] > -INF else (c["ub"][j] - 6 if c["ub"][j] < INF else -5.0)
        hi = c["ub"][j] if c["ub"][j] < INF else lo + 6
        if c["is_int"][j]:
            x0.append(float(rng.randint(int(math.ceil(lo)), int(math.floor(hi)))))
        else:
            x0.append(lo + 0.5 * rng.randint(0, int(2 * (hi - lo))))
    coefs = [1, 1, 1, 2, 3, -1, -1, -2, 0.5, 1.5, 4, 7]
    for i in range(m + (1 if rng.random() < 0.8 else 0)):
        if i == m:
            # Bounding row: plain sum in a 2-sided range keeps most cases
            # bounded while staying random.
            c["rows"].append([(j, 1.0) for j in range(n)])
            c["senses"].append("L")
            c["rhs"].append(0.0)
            c["ranges"].append(0.0)
            c["m"] += 1
        else:
            c["rows"][i] = [(j, float(rng.choice(coefs))) for j in range(n)
                            if rng.random() < 0.55] or [(rng.randrange(n), 1.0)]
            c["senses"][i] = rng.choices("LGE", (0.55, 0.3, 0.15))[0]
        if not planted:
            c["rhs"][i] = float(rng.choice([0, 2, 3, 5, 7, 10, 12, -3, 4.5, 20]))
            if rng.random() < 0.12:
                c["ranges"][i] = float(rng.choice([1.5, 2.0, 4.0, -2.0, -3.0, 6.0]))
            continue
        act = sum(a * x0[j] for j, a in c["rows"][i])
        slack = 0.0 if c["senses"][i] == "E" else float(rng.choice([0, 0, 1, 2.5, 5, 9]))
        if i == m:
            slack = float(rng.randint(0, 10))
        c["rhs"][i] = act + slack if c["senses"][i] == "L" else act - slack
        if i == m or rng.random() < 0.12:
            w = slack + float(rng.choice([0, 1.5, 3, 20 if i == m else 6]))
            c["ranges"][i] = w if c["senses"][i] != "E" else rng.choice([w, -w])
    c["planted"] = x0 if planted else None
    return c


def _parity_infeasible(rng, k, lp_unbounded):
    """sum_j 2 a_j x_j = odd with all-integer x: no integer point (parity).

    lp_unbounded=False: all x bounded, LP relaxation feasible and bounded.
    lp_unbounded=True:  x free integers with an objective along the
    equality's null space -> LP relaxation unbounded, MIP infeasible. (With
    rational data that is the only way an LP-unbounded MIP can be bounded:
    Meyer 1974, "On the existence of optimal solutions to integer and mixed-
    integer programming problems", Math. Prog. 7.)
    """
    n = rng.randint(3, 8)
    m = rng.randint(1, 4)
    c = _empty_case(f"rand_{k:03d}_parity{'_lpunb' if lp_unbounded else ''}", n, m + 1)
    c["maximize"] = rng.random() < 0.5
    a = [2 * rng.choice([1, -1, 2, 3]) for _ in range(n)]
    c["rows"][0] = [(j, float(a[j])) for j in range(n)]
    c["senses"][0] = "E"
    c["rhs"][0] = float(2 * rng.randint(-3, 6) + 1)
    for i in range(1, m + 1):
        c["rows"][i] = [(j, float(rng.choice([1, 2, -1, 3]))) for j in range(n)
                        if rng.random() < 0.5] or [(0, 1.0)]
        c["senses"][i] = "L"
        c["rhs"][i] = float(rng.randint(10, 30))
        if lp_unbounded:
            c["rows"][i] = [(j, 0.0) for j in range(n)]  # placeholder, dropped below
    for j in range(n):
        c["is_int"][j] = True
        if lp_unbounded:
            c["lb"][j], c["ub"][j] = -INF, INF
        else:
            c["lb"][j], c["ub"][j] = float(rng.choice([0, -3])), float(rng.choice([4, 6, 9]))
        c["obj"][j] = float(rng.choice([-2, -1, 1, 2, 3]))
    if lp_unbounded:
        # Only the parity row: x0 and x1 move along a0*d0 + a1*d1 = 0.
        c["m"] = 1
        c["rows"], c["senses"], c["rhs"], c["ranges"] = (c["rows"][:1], c["senses"][:1],
                                                         c["rhs"][:1], [0.0])
        # objective strictly improves along d = (a1, -a0, 0, ...)
        sign = 1.0 if c["maximize"] else -1.0
        c["obj"] = [0.0] * n
        c["obj"][0] = sign * (1.0 if a[1] > 0 else -1.0)
    c["expect"] = "Infeasible"
    c["family"] = "random_parity"
    return c


def _unbounded_mip(rng, k):
    """Feasible MIP with an integer recession direction -> MIP unbounded."""
    n = rng.randint(3, 8)
    m = rng.randint(1, 5)
    c = _empty_case(f"rand_{k:03d}_unbounded", n, m)
    c["maximize"] = rng.random() < 0.5
    for j in range(n):
        c["is_int"][j] = j == 0 or rng.random() < 0.6
        c["lb"][j], c["ub"][j] = 0.0, float(rng.choice([1, 3, 10]))
        c["obj"][j] = float(rng.choice([-1, 1, 2]))
    # x0 >= 0 unbounded above, appears in rows only with signs that never
    # block growth; objective improves with x0.
    c["ub"][0] = INF
    c["obj"][0] = 1.0 if c["maximize"] else -1.0
    for i in range(m):
        c["rows"][i] = [(j, float(rng.choice([1, 2, 3]))) for j in range(1, n)
                        if rng.random() < 0.6] or [(1, 1.0)]
        c["rows"][i].append((0, -1.0))  # x0 only relaxes a <= row
        c["senses"][i] = "L"
        c["rhs"][i] = float(rng.randint(2, 15))
    c["expect"] = "Unbounded"
    c["family"] = "random_unbounded"
    return c


def model_to_case(mdl: "G.Model", fam: str) -> dict:
    n, m = mdl.ncols, mdl.nrows
    c = _empty_case(mdl.name, n, m)
    c["family"] = fam
    c["maximize"] = mdl.obj_sense == "maximize"
    c["bv"] = True
    for j, v in enumerate(mdl.vars):
        c["lb"][j], c["ub"][j] = float(v.lb), float(v.ub)
        c["is_int"][j] = v.vtype != G.VTYPE_CONT
    for j, a in mdl.obj:
        c["obj"][int(j)] += float(a)
    for i, r in enumerate(mdl.rows):
        c["rows"][i] = [(int(j), float(a)) for j, a in r.coeffs]
        c["senses"][i] = {"<=": "L", ">=": "G", "==": "E"}[r.sense]
        c["rhs"][i] = float(r.rhs)
    p = mdl.meta.get("planted")
    if p is not None:
        x = [0.0] * n
        for j, v in (p.items() if isinstance(p, dict) else enumerate(p)):
            x[int(j)] = float(v)
        c["planted"] = x
    if mdl.meta.get("infeasible"):
        c["expect"] = "Infeasible"
    return c


# family label -> extra kwargs on top of MIP_SMALL
VARIANTS = [(f, f, {}) for f in G.MIP_FAMILIES] + [
    ("uncap_facility_location[weak]", "uncap_facility_location", {"form": "weak"}),
    ("gen_int_knapsack[eq]", "gen_int_knapsack", {"equality": True}),
    ("gen_int_knapsack[infeasible]", "gen_int_knapsack", {"infeasible": True}),
    ("set_partition[infeasible]", "set_partition", {"infeasible": True}),
    ("gap[infeasible]", "gap", {"infeasible": True}),
    # wide-coefficient / objective-shape / mixed-structure variants
    ("lot_sizing[bigM]", "lot_sizing", {"big_m": 1e5}),
    ("fixed_charge_network[bigM]", "fixed_charge_network", {"big_m": 1e4}),
    ("fixed_charge_network[fixed]", "fixed_charge_network", {"var_cost": False}),
    ("jobshop_bigM[wideM]", "jobshop_bigM", {"big_m": 1e6, "p_max": 100}),
    ("network_design[fixed]", "network_design", {"var_cost": False}),
    ("network_design[modular]", "network_design", {"modules": 3}),
    ("int_production[narrow]", "int_production", {"wide": False}),
    ("block_angular[nolink]", "block_angular", {"link_rows": 0}),
    ("block_angular[infeasible]", "block_angular", {"infeasible": True}),
    ("general_ip[min]", "general_ip", {"maximize": False}),
    ("lights_out[infeasible]", "lights_out", {"infeasible": True}),
    ("p_median[bounded]", "p_median", {"bounded": True}),
]


def family_cases(seeds, fam_filter=None):
    out = []
    for label, fam, kw in VARIANTS:
        if fam_filter and fam_filter not in (label, fam):
            continue
        for s in seeds:
            mdl = G.REGISTRY[fam](**{**G.MIP_SMALL[fam], **kw, "seed": s})
            err = G.validate_planted(mdl, eps=1e-6)
            if err:
                raise SystemExit(f"GENERATOR ASSERTION FAILED {label} seed {s}: {err}")
            c = model_to_case(mdl, label)
            c["name"] = f"{label.replace('[', '_').replace(']', '')}_s{s}"
            out.append(c)
    return out


# --------------------------------------------------------------------------- oracles


def solve_highs(mdl, time_limit=ORACLE_TIME, presolve=True):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("mip_rel_gap", 0.0)
    h.setOptionValue("mip_abs_gap", 1e-9)
    h.setOptionValue("threads", 1)
    h.setOptionValue("random_seed", 0)
    h.setOptionValue("time_limit", time_limit)
    if not presolve:
        h.setOptionValue("presolve", "off")
    n = mdl.n
    h.addVars(n, list(mdl.lb), list(mdl.ub))
    if n:
        h.changeColsCost(n, list(range(n)), list(mdl.obj))
        h.changeColsIntegrality(n, list(range(n)),
                                [highspy.HighsVarType.kInteger if t else
                                 highspy.HighsVarType.kContinuous for t in mdl.is_int])
    h.changeObjectiveOffset(mdl.obj_offset)
    h.changeObjectiveSense(highspy.ObjSense.kMaximize if mdl.maximize
                           else highspy.ObjSense.kMinimize)
    starts, index, value = [], [], []
    for row in mdl.rows:
        starts.append(len(index))
        for j, a in row:
            index.append(j)
            value.append(a)
    if mdl.m:
        h.addRows(mdl.m, list(mdl.row_lo), list(mdl.row_hi), len(index), starts, index, value)
    h.run()
    st = h.getModelStatus()
    status = {HMS.kOptimal: "Optimal", HMS.kInfeasible: "Infeasible",
              HMS.kUnbounded: "Unbounded", HMS.kUnboundedOrInfeasible: IORU,
              HMS.kTimeLimit: "TimeLimit"}.get(st, f"Other({st.name})")
    obj, x = None, None
    if status == "Optimal":
        obj = h.getInfo().objective_function_value
        x = list(h.getSolution().col_value)
    return status, obj, x


def solve_scip(mdl, time_limit=ORACLE_TIME):
    s = pyscipopt.Model()
    s.hideOutput()
    s.setParam("limits/gap", 0.0)
    s.setParam("limits/absgap", 1e-9)
    s.setParam("limits/time", time_limit)
    s.setParam("randomization/randomseedshift", 0)
    xs = []
    for j in range(mdl.n):
        lo = None if mdl.lb[j] == -INF else mdl.lb[j]
        hi = None if mdl.ub[j] == INF else mdl.ub[j]
        xs.append(s.addVar(name=f"x{j}", vtype="I" if mdl.is_int[j] else "C", lb=lo, ub=hi,
                           obj=mdl.obj[j]))
    for i, row in enumerate(mdl.rows):
        lo, hi = mdl.row_lo[i], mdl.row_hi[i]
        if not row:
            if lo > 1e-9 or hi < -1e-9:  # empty row 0 outside [lo, hi]
                s.addCons(pyscipopt.scip.ExprCons(pyscipopt.Expr() + 0.0 * xs[0], lhs=1.0))
            continue
        expr = pyscipopt.quicksum(a * xs[j] for j, a in row)
        s.addCons(pyscipopt.scip.ExprCons(expr, lhs=None if lo == -INF else lo,
                                          rhs=None if hi == INF else hi))
    s.addObjoffset(mdl.obj_offset)
    if mdl.maximize:
        s.setMaximize()
    s.optimize()
    st = s.getStatus()
    status = {"optimal": "Optimal", "infeasible": "Infeasible", "unbounded": "Unbounded",
              "inforunbd": IORU, "timelimit": "TimeLimit"}.get(st, f"Other({st})")
    obj, x = None, None
    if status == "Optimal":
        obj = s.getObjVal()
        sol = s.getBestSol()
        x = [sol[v] for v in xs]
    return status, obj, x


def solve_scip_guarded(mdl, timeout=ORACLE_TIME + 30):
    """solve_scip in a forked child with a hard kill.

    SCIP 10 can ignore limits/time: with numerics/feastol=1e-9 it spun
    forever on an unbounded random case (seed 2, rand_111), so the child is
    killed after `timeout` and the case reported as 'Hang'.
    """
    ctx = multiprocessing.get_context("fork")
    rd, wr = ctx.Pipe(duplex=False)
    proc = ctx.Process(target=lambda: wr.send(solve_scip(mdl)))
    proc.start()
    res = rd.recv() if rd.poll(timeout) else ("Hang", None, None)
    proc.kill()
    proc.join()
    return res


def obj_close(a, b, tol=1e-6):
    return a is not None and b is not None and abs(a - b) <= tol * max(1.0, abs(b))


# SCIP runs at its default feastol 1e-6 and legitimately drifts continuous
# values by ~1e-6, moving the optimum by ~1e-6 relative (rand_173/285 at seed
# 20261022); the ORACLE cross-check therefore uses 1e-5. Our solver is still
# held to 1e-6 against HiGHS.
SCIP_OBJ_TOL = 1e-5


def oracles_agree(hs, hobj, ss, sobj):
    if ss == "Hang":
        return True  # no cross-check available; counted separately
    if hs == ss == "Optimal":
        return obj_close(sobj, hobj, SCIP_OBJ_TOL)
    term = {"Infeasible", "Unbounded", IORU}
    if hs in term and ss in term:
        return hs == ss or IORU in (hs, ss)
    return False


# --------------------------------------------------------------------------- our solver


def run_ours(binary, path, timeout, extra):
    cmd = [binary, "mip", path, "--csv"] + extra
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return dict(status="Timeout", obj=None, x=None, wall=time.time() - t0, err="timeout")
    wall = time.time() - t0
    lines = [ln for ln in p.stdout.splitlines() if ln.strip()]
    if len(lines) < 2 or not lines[0].startswith("status"):
        return dict(status="ERR", obj=None, x=None, wall=wall,
                    err=(p.stderr.strip().splitlines() or ["no csv"])[-1][:80])
    parts = lines[1].split(",", 6)
    try:
        obj = float(parts[1])
    except ValueError:
        obj = math.nan
    try:
        x = parse_vector(parts[6]) if len(parts) > 6 else None
    except ValueError:
        x = None
    return dict(status=parts[0], obj=obj, x=x, wall=wall, err="",
                nodes=parts[3] if len(parts) > 3 else "", exit=p.returncode)


# --------------------------------------------------------------------------- io self-check


def roundtrip_problems(case, mdl):
    """Our MPS writer + independent reader must reproduce the case exactly."""
    if (mdl.n, mdl.m) != (case["n"], case["m"]):
        return f"shape {mdl.n}x{mdl.m} != {case['n']}x{case['m']}"
    for j in range(mdl.n):
        if (mdl.lb[j], mdl.ub[j], mdl.is_int[j], mdl.obj[j]) != (
                case["lb"][j], case["ub"][j], case["is_int"][j], case["obj"][j]):
            return f"col {j} mismatch"
    for i in range(mdl.m):
        want = sorted((j, a) for j, a in case["rows"][i] if a != 0.0)
        if sorted(mdl.rows[i]) != want:
            return f"row {i} coeffs mismatch"
        if (mdl.row_lo[i], mdl.row_hi[i]) != row_bounds(case["senses"][i], case["rhs"][i],
                                                       case["ranges"][i]):
            return f"row {i} bounds mismatch"
    if mdl.maximize != case["maximize"]:
        return "sense mismatch"
    return None


def rusolve_dump_problems(binary, path, mdl):
    """rusolve-io must read the same bounds/rows (integrality is not dumped)."""
    p = subprocess.run([binary, "dump", path], capture_output=True, text=True)
    if p.returncode != 0:
        return f"dump failed: {p.stderr.strip()[:80]}"
    nz, bnd = {}, {}
    for ln in p.stdout.splitlines():
        t = ln.split()
        if t[0] == "nz":
            nz[(int(t[2]), int(t[1]))] = float(t[3])
        elif t[0] == "bound":
            bnd[int(t[1])] = (float(t[2]), float(t[3]))
    want = {(i, j): a for i, row in enumerate(mdl.rows) for j, a in row if a != 0.0}
    if nz != want:
        return "dump nz mismatch"
    for j in range(mdl.n):
        if bnd.get(j) != (mdl.lb[j], mdl.ub[j]):
            return f"dump bound {j} {bnd.get(j)} != {(mdl.lb[j], mdl.ub[j])}"
    return None


# --------------------------------------------------------------------------- main


def fmt(v):
    return "-" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.6g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", type=int, default=80, help="number of random cases (a)")
    ap.add_argument("--seed", type=int, default=20261022)
    ap.add_argument("--family-seeds", type=int, default=4, help="seeds per family variant (b)")
    ap.add_argument("--family", default=None, help="only this family / variant label")
    ap.add_argument("--only-random", action="store_true")
    ap.add_argument("--oracle-only", action="store_true",
                    help="validate generators + checker + oracles without our binary")
    ap.add_argument("--timeout", type=float, default=30.0, help="per-case limit for our binary")
    ap.add_argument("--extra", default="", help="extra args for `rusolve-cli mip`")
    ap.add_argument("--cases-dir", default=None,
                    help="where case .mps files go (default: fresh run-<time>-<pid> dir "
                         "under tools/diff/mip_cases/, so concurrent gates never clobber)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    binary = os.environ.get("RUSOLVE_BIN", os.path.join(REPO, "target", "release", "rusolve-cli"))
    have_bin = os.path.exists(binary)
    if not args.oracle_only and not have_bin:
        raise SystemExit(f"{binary} missing (build it, or pass --oracle-only)")

    cases = []
    if args.family is None or args.family.startswith("random"):
        rng = random.Random(args.seed)
        cases += [gen_random_case(rng, k) for k in range(args.cases)]
        if args.family:
            cases = [c for c in cases if c["family"] == args.family]
    if not args.only_random and not (args.family or "").startswith("random"):
        seeds = [args.seed % 1000 + s for s in range(args.family_seeds)]
        cases += family_cases(seeds, args.family)
    if not cases:
        raise SystemExit("no cases selected")

    if args.cases_dir:
        cases_dir = args.cases_dir
        os.makedirs(cases_dir, exist_ok=True)
    else:
        os.makedirs(CASES, exist_ok=True)
        cases_dir = tempfile.mkdtemp(prefix=time.strftime("run-%Y%m%d-%H%M%S-")
                                     + f"{os.getpid()}-", dir=CASES)
    rows, n_pass, n_fail, n_dis, n_hang, n_tie, n_nonterm = [], 0, 0, 0, 0, 0, 0
    for c in cases:
        path = os.path.join(cases_dir, c["name"] + ".mps")
        write_mps(c, path)
        mdl = read_mps(path)
        notes = []
        p = roundtrip_problems(c, mdl)
        if p:
            notes.append("ROUNDTRIP " + p)
        hs, hobj, hx = solve_highs(mdl)
        ss, sobj, sx = solve_scip_guarded(mdl)
        disagree = not oracles_agree(hs, hobj, ss, sobj)
        n_hang += ss == "Hang"
        info = []
        if disagree:
            # 2-of-3 tiebreak: HiGHS 1.15.1 presolve returns wrong statuses on
            # rare random cases (seed 3: rand_027 Optimal on an unbounded MIP,
            # rand_318 Infeasible with a valid planted point); presolve-off
            # HiGHS + SCIP agreeing becomes the reference, flagged in the table.
            ps, pobj, px = solve_highs(mdl, presolve=False)
            if oracles_agree(ps, pobj, ss, sobj):
                info.append(f"HIGHS-PRESOLVE-BUG (presolve on: {hs} {fmt(hobj)}; "
                            f"reference = presolve off + SCIP)")
                hs, hobj, hx = ps, pobj, px
                disagree = False
                n_tie += 1
        if hx is not None:
            r = check_solution(mdl, hx, hobj)
            if not r.ok:
                notes.append("HIGHS-X " + r.short())
        if c["expect"] and hs not in (c["expect"], IORU):
            notes.append(f"EXPECT {c['expect']}")
        verdict = "PASS"
        ours = None
        if args.oracle_only:
            if sx is not None:
                r = check_solution(mdl, sx, sobj)
                if not r.ok:
                    notes.append("SCIP-X " + r.short())
            if c["planted"] is not None:
                r = check_solution(mdl, c["planted"])
                if not r.ok:
                    notes.append("PLANTED " + r.short())
                elif hobj is not None and (r.obj - hobj) * (1 if mdl.maximize else -1) > \
                        1e-6 * max(1.0, abs(hobj)):
                    notes.append(f"PLANTED beats HiGHS opt ({r.obj} vs {hobj})")
                elif hs == "Infeasible":
                    notes.append("PLANTED feasible but HiGHS Infeasible")
            if have_bin:
                p = rusolve_dump_problems(binary, path, mdl)
                if p:
                    notes.append("RUSOLVE-IO " + p)
            if notes:
                verdict = "FAIL"
        else:
            # Our own --time-limit sits below the hard kill so a clean
            # TimeLimit row (with bound/nodes) is the normal timeout path.
            tl = ["--time-limit", f"{max(1.0, 0.8 * args.timeout):g}"]
            ours = run_ours(binary, path, args.timeout, tl + args.extra.split())
            st = ours["status"]
            if st in ("Timeout", "TimeLimit", "NodeLimit") and hs in ("Infeasible", IORU) \
                    and any(mdl.is_int[j] and (mdl.lb[j] == -INF or mdl.ub[j] == INF)
                            for j in range(mdl.n)):
                # Known limitation: plain B&B cannot prove infeasibility over
                # unbounded integer boxes (needs gcd reasoning). Reported, not
                # failed.
                verdict = "EXPECTED-NONTERMINATION"
            elif st in ("Timeout", "TimeLimit", "NodeLimit"):
                verdict = f"FAIL({st.lower()})"
            elif hs == "Optimal":
                if st != "Optimal":
                    verdict = "FAIL(status)"
                elif not obj_close(ours["obj"], hobj):
                    verdict = "FAIL(obj)"
                if st == "Optimal":
                    if ours["x"] is None:
                        verdict = "FAIL(no-x)"
                    else:
                        r = check_solution(mdl, ours["x"], ours["obj"])
                        if not r.ok:
                            verdict = "FAIL(check)" if verdict == "PASS" else verdict
                            notes.append("OUR-X " + r.short())
            elif hs == IORU:
                if st not in ("Infeasible", "Unbounded"):
                    verdict = "FAIL(status)"
            elif st != hs:
                verdict = "FAIL(status)"
                if st == "Optimal" and ours["x"] is not None:
                    r = check_solution(mdl, ours["x"], ours["obj"])
                    notes.append("our x " + ("FEASIBLE (oracle suspect)" if r.ok else r.short()))
            if ours.get("err"):
                notes.append(ours["err"])
            if any(nt.startswith(("ROUNDTRIP", "HIGHS-X", "EXPECT ")) for nt in notes):
                verdict = "FAIL(harness)"
        if disagree:
            verdict = "ORACLE-DISAGREE"
            notes.append(f"scip {ss} {fmt(sobj)}")
        if verdict == "PASS":
            n_pass += 1
        elif verdict == "ORACLE-DISAGREE":
            n_dis += 1
        elif verdict == "EXPECTED-NONTERMINATION":
            n_nonterm += 1
        else:
            n_fail += 1
        rows.append((c, hs, hobj, ss, ours, verdict, notes + info))

    w = max(len(c["name"]) for c in cases)
    if args.oracle_only:
        print(f"{'case':<{w}} {'n':>3} {'m':>3} {'highs':>11} {'obj':>12} {'scip':>11}  verdict")
    else:
        print(f"{'case':<{w}} {'n':>3} {'m':>3} {'ours':>11} {'obj':>12} {'highs':>11} "
              f"{'obj':>12} {'t(s)':>6}  verdict")
    for c, hs, hobj, ss, ours, verdict, notes in rows:
        if not args.verbose and verdict == "PASS" and len(rows) > 60 and not notes \
                and ss != "Hang":
            continue
        head = f"{c['name']:<{w}} {c['n']:>3} {c['m']:>3} "
        if args.oracle_only:
            line = head + f"{hs:>11} {fmt(hobj):>12} {ss:>11}  {verdict}"
        else:
            line = head + (f"{ours['status']:>11} {fmt(ours['obj']):>12} {hs:>11} "
                           f"{fmt(hobj):>12} {ours['wall']:>6.2f}  {verdict}")
        if notes:
            line += "  " + "; ".join(notes)
        print(line)

    by_fam = {}
    for c, hs, *_rest in rows:
        f = by_fam.setdefault(c["family"], {})
        f[hs] = f.get(hs, 0) + 1
    print("\ncases per family (HiGHS status counts):")
    for f, d in by_fam.items():
        print(f"  {f:<34} {sum(d.values()):>4}  " + " ".join(f"{k}={v}" for k, v in sorted(d.items())))
    mode = "oracle-only" if args.oracle_only else os.path.relpath(binary, REPO)
    print(f"\n[{mode}] {len(rows)} cases: {n_pass} PASS, {n_fail} FAIL, {n_dis} ORACLE-DISAGREE"
          + (f", {n_nonterm} EXPECTED-NONTERMINATION" if n_nonterm else "")
          + (f" ({n_hang} without SCIP cross-check: SCIP hang)" if n_hang else "")
          + (f" ({n_tie} HiGHS-presolve disagreements resolved 2-of-3)" if n_tie else ""))
    print(f"case files: {os.path.relpath(cases_dir, REPO)}")
    all_pass = n_fail == 0 and n_dis == 0
    print("ALL PASS" if all_pass else "FAILURES PRESENT")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
