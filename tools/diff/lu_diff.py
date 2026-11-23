#!/usr/bin/env python3
"""Factorization differential gate: rusolve LP objective vs highspy (rel 1e-6).

Solves every `.mps` in the given directories (default: benches/data/netlib)
with `target/release/rusolve-cli solve --csv` (primal and dual, each with presolve on AND `--no-presolve`:
branch-and-bound solves presolve-off LPs, and a presolve-on-only check once missed a
silent primal wrong answer on tuff) and
compares the objective to HiGHS simplex. Exit 1 on any objective mismatch
or non-optimal status where HiGHS is optimal. Known iteration-limit modes
live in KNOWN_MODES (empty at present).

Usage: python tools/diff/lu_diff.py [dir ...] [--only a,b]
Covers the sparse LU factorization and its refactorization trigger.
"""
import os
import subprocess
import sys
import time

import highspy

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CLI = os.path.join(REPO, "target", "release", "rusolve-cli")
KNOWN = set()
# The former primal iteration-limit modes (brandy primal without presolve,
# e226 primal with presolve) are fixed; an entry reappearing here is a
# REGRESSION, not a known limit.
KNOWN_MODES = set()  # (name, dual, presolve)


def highs_obj(path):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("solver", "simplex")
    h.readModel(path)
    h.run()
    ok = h.getModelStatus() == highspy.HighsModelStatus.kOptimal
    return ok, h.getInfo().objective_function_value


def ours(path, dual, presolve=True):
    env = dict(os.environ)
    env.pop("RUSOLVE_PRIMAL", None)
    if not dual:
        env["RUSOLVE_PRIMAL"] = "1"  # dual is the CLI default (dual-default-eval §9)
    cmd = [CLI, "solve", path, "--csv"] + ([] if presolve else ["--no-presolve"])
    p = subprocess.run(cmd, capture_output=True, text=True,
                       timeout=300, env=env)
    ln = [l for l in p.stdout.splitlines() if l and not l.startswith("status")]
    if not ln:
        return "NOOUTPUT", float("nan"), 0
    s, obj, p1, p2 = ln[-1].split(",")[:4]
    return s, float(obj), int(p1) + int(p2)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    only = None
    if "--only" in sys.argv:
        only = set(sys.argv[sys.argv.index("--only") + 1].split(","))
        args = [a for a in args if "," not in a and a not in only]
    dirs = args or [os.path.join(REPO, "benches", "data", "netlib")]
    bad = 0
    for d in dirs:
        for f in sorted(os.listdir(d)):
            if not f.endswith(".mps"):
                continue
            name = f[:-4]
            if only and name not in only:
                continue
            path = os.path.join(d, f)
            ok, ho = highs_obj(path)
            if not ok:
                continue
            for dual, presolve in ((False, True), (True, True), (False, False), (True, False)):
                t = time.perf_counter()
                s, o, it = ours(path, dual, presolve)
                dt = time.perf_counter() - t
                good = s == "Optimal" and abs(o - ho) <= 1e-6 * max(1.0, abs(ho))
                known = name in KNOWN or (name, dual, presolve) in KNOWN_MODES
                tag = "ok" if good else ("KNOWN" if known else "FAIL")
                if not good and not known:
                    bad += 1
                print(f"{name:24s} {'dual' if dual else 'prim'}{'' if presolve else '-np'} {s:10s} {o:.8g} vs {ho:.8g} "
                      f"it={it} {dt:.2f}s {tag}", flush=True)
    print("FAILURES:", bad)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
