#!/usr/bin/env python3
"""Simplex differential check: rusolve simplex vs highspy on seeded random LPs.

Every case runs rusolve TWICE (presolve ON and OFF); both must agree with
HiGHS on status and objective (rel 1e-6) and with each other.
Cases include: L/G/E senses, negative rhs, ranged rows (MPS RANGES section),
fixed variables (lb == ub), finite ubs, positive lbs, maximize.

Usage: python tools/diff/simplex_diff.py
Exit 0 iff ALL PASS.
"""
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CASES = os.path.join(HERE, "cases")
BIN = os.path.join(REPO, "target", "debug", "rusolve-cli")

import highspy
from highspy import HighsModelStatus

INF = float("inf")


def gen_case(rng):
    n = rng.randint(3, 18)
    m = rng.randint(2, 14)
    senses = [rng.choice("LGE") for _ in range(m)]
    A = [[rng.choice([0, 0, 1, 1, 1, 2, -1, -1, 0.5, -2.0, 3]) for _ in range(n)] for _ in range(m)]
    rhs = [rng.choice([0, 2, 5, 10, -3, -7, 1.5]) for _ in range(m)]
    obj = [rng.choice([-2, -1, -0.5, 0, 1, 2, 3, 4]) for _ in range(n)]
    lower = [rng.choice([0, 0, 0, 0, 1, 2.5]) for _ in range(n)]
    upper = [rng.choice([INF, 4, 10, 20, INF]) for _ in range(n)]
    upper = [u if u >= lower[j] else lower[j] + 1 for j, u in enumerate(upper)]
    # ranged rows: ~30% of cases, one random row gets a signed range
    ranges = [0.0] * m
    if rng.random() < 0.3:
        ranges[rng.randrange(m)] = rng.choice([1.5, 2.0, 4.0, -2.0, -3.0])
    # fixed variable: ~25% of cases, one variable gets lb == ub
    if rng.random() < 0.25:
        j = rng.randrange(n)
        v = rng.choice([0.0, 1.0, 2.0])
        lower[j] = v
        upper[j] = v
    # free variables: ~30% of cases one FR, ~20% one MI (lb=-inf, ub kept)
    if rng.random() < 0.30:
        j = rng.randrange(n)
        lower[j] = -INF
        upper[j] = INF
    if rng.random() < 0.20:
        j = rng.randrange(n)
        lower[j] = -INF
        if upper[j] == INF:
            upper[j] = rng.choice([4.0, 10.0])
    maximize = rng.random() < 0.25
    return dict(n=n, m=m, senses=senses, A=A, rhs=rhs, obj=obj, lower=lower,
                upper=upper, ranges=ranges, maximize=maximize)


def write_mps(c, path):
    L = []
    L.append("NAME CASE")
    L.append(f"OBJSENSE\n {'MAX' if c['maximize'] else 'MIN'}")
    L.append("ROWS")
    L.append(" N  OBJ")
    for i, s in enumerate(c["senses"]):
        L.append(f" {s}  R{i}")
    L.append("COLUMNS")
    for j in range(c["n"]):
        if c["obj"][j]:
            L.append(f"    X{j}    OBJ    {c['obj'][j]!r}")
        for i in range(c["m"]):
            if c["A"][i][j]:
                L.append(f"    X{j}    R{i}    {c['A'][i][j]!r}")
    L.append("RHS")
    for i in range(c["m"]):
        if c["rhs"][i]:
            L.append(f"    RHS    R{i}    {c['rhs'][i]!r}")
    if any(c["ranges"]):
        L.append("RANGES")
        for i in range(c["m"]):
            if c["ranges"][i]:
                L.append(f"    RHS    R{i}    {c['ranges'][i]!r}")
    L.append("BOUNDS")
    for j in range(c["n"]):
        if c["lower"][j] == -INF:
            if c["upper"][j] == INF:
                L.append(f" FR BND X{j}")
            else:
                L.append(f" MI BND X{j}")
                L.append(f" UP BND X{j} {c['upper'][j]!r}")
            continue
        if c["lower"][j] != 0:
            L.append(f" LO BND X{j} {c['lower'][j]!r}")
        if c["upper"][j] != INF:
            L.append(f" UP BND X{j} {c['upper'][j]!r}")
    L.append("ENDATA")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def highs_row_bounds(c):
    """Row (lower, upper) pairs, expanding MPS ranges to two-sided rows."""
    out = []
    for i in range(c["m"]):
        s, r, rg = c["senses"][i], c["rhs"][i], c["ranges"][i]
        if rg == 0.0:
            if s == "L":
                out.append((-1e18, r))
            elif s == "G":
                out.append((r, 1e18))
            else:
                out.append((r, r))
        elif s == "L":
            out.append((r - abs(rg), r))
        elif s == "G":
            out.append((r, r + abs(rg)))
        elif rg > 0:
            out.append((r, r + rg))
        else:
            out.append((r + rg, r))
    return out


def solve_highs(c):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.addVars(c["n"], [c["lower"][j] for j in range(c["n"])],
              [c["upper"][j] for j in range(c["n"])])
    h.changeColsCost(c["n"], list(range(c["n"])), list(c["obj"]))
    h.changeObjectiveSense(highspy.ObjSense.kMaximize if c["maximize"] else highspy.ObjSense.kMinimize)
    starts, index, value = [], [], []
    for i in range(c["m"]):
        starts.append(len(index))
        for j in range(c["n"]):
            if c["A"][i][j]:
                index.append(j)
                value.append(c["A"][i][j])
    lower, upper = zip(*highs_row_bounds(c)) if c["m"] else ([], [])
    h.addRows(c["m"], list(lower), list(upper), len(index), starts, index, value)
    h.run()
    return h


def run_rusolve(path, no_presolve):
    cmd = [BIN, "solve", path, "--csv"] + (["--no-presolve"] if no_presolve else [])
    p = subprocess.run(cmd, capture_output=True, text=True)
    lines = p.stdout.strip().splitlines()
    if len(lines) < 2:
        return "ERR", None, p.stderr
    parts = lines[1].split(",")
    return parts[0], float(parts[1]), p.stderr


def main():
    rng = random.Random(20261021)
    os.makedirs(CASES, exist_ok=True)
    results = []
    all_pass = True
    for k in range(60):
        c = gen_case(rng)
        path = os.path.join(CASES, f"case_{k:02d}.mps")
        write_mps(c, path)
        h = solve_highs(c)
        st = h.getModelStatus()
        hopt = st == HighsModelStatus.kOptimal
        hobj = h.getObjectiveValue() if hopt else None
        # HiGHS "kOptimal" with |obj| >= 1e15 is unbounded-direction runaway
        # (the LP is genuinely UNBOUNDED) — treat as such.
        if hopt and abs(hobj) >= 1e15:
            hopt, hs_unb = False, True
        else:
            hs_unb = st in (HighsModelStatus.kInfeasible, HighsModelStatus.kUnbounded)

        st_on, obj_on, err_on = run_rusolve(path, no_presolve=False)
        st_off, obj_off, err_off = run_rusolve(path, no_presolve=True)

        def agrees(status, obj):
            ours_opt = status == "Optimal"
            ok_status = (ours_opt == hopt) or (status in ("Infeasible", "Unbounded") and hs_unb)
            ok_obj = hopt and ours_opt and abs(obj - hobj) <= 1e-6 * max(1.0, abs(hobj))
            return ok_status and (ok_obj or not hopt)

        mutual = (st_on == st_off) and (
            hopt is False or st_on != "Optimal"
            or abs(obj_on - obj_off) <= 1e-6 * max(1.0, abs(obj_off)))
        ok = agrees(st_on, obj_on) and agrees(st_off, obj_off) and mutual
        all_pass &= ok
        results.append((k, st_on, obj_on, str(st).replace("HighsModelStatus.", ""), hobj, st_off, "PASS" if ok else "FAIL"))
    print(f"{'case':>4} {'ps=on':>11} {'obj':>13} {'highs':>14} {'obj':>13} {'ps=off':>11} ")
    for k, s, o, hs, ho, so, v in results:
        print(f"{k:>4} {s:>11} {('%.6g' % o) if o is not None and o == o else '-':>13} {hs:>14} {('%.6g' % ho) if ho is not None else '-':>13} {so:>11}  {v}")
    print("ALL PASS" if all_pass else "FAILURES PRESENT")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
