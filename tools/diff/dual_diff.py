#!/usr/bin/env python3
"""Dual-simplex differential gate: legacy dual (RUSOLVE_EKK=0) objective vs highspy (rel 1e-6).

Covers the HiGHS-derived pieces of the dual engine (up-front cost
perturbation, dual steepest edge, plateau test, gated Harris): Netlib with presolve on and off, plus
a handful of the 47 relaxed MIP-root LPs when a directory of them is given.
Also reports the "primal fallback" count (silent dual bails).

Usage: python tools/diff/dual_diff.py [--relaxed DIR]
Exit 1 on any wrong status/objective where HiGHS is optimal.
"""
import os
import subprocess
import sys

import highspy

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CLI = os.environ.get("RUSOLVE_BIN", os.path.join(REPO, "target", "release", "rusolve-cli"))
RELAXED = ["graph20-20-1rand", "neos-3004026-krka", "milo-v12-6-r2-40-1",
           "drayage-25-23", "ns1830653"]
KNOWN = set()  # instances with a known mismatch (none at present)


def highs_obj(path):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("solver", "simplex")
    h.readModel(path)
    h.run()
    return h.getModelStatus() == highspy.HighsModelStatus.kOptimal, h.getInfo().objective_function_value


def ours(path, extra):
    # RUSOLVE_EKK=0 selects the legacy engine, whose cold default is the dual; the old
    # RUSOLVE_DUAL=1 is not a recognised variable (the CLI rejects it)
    env = dict(os.environ, RUSOLVE_EKK="0")
    p = subprocess.run([CLI, "solve", path, "--csv"] + extra, capture_output=True, text=True,
                       timeout=300, env=env)
    ln = [l for l in p.stdout.splitlines() if l and not l.startswith("status")]
    s, obj, p1, p2 = ln[-1].split(",")[:4]
    return s, float(obj), int(p1) + int(p2), "primal fallback" in p.stderr


def main():
    jobs = []
    nl = os.path.join(REPO, "benches", "data", "netlib")
    for f in sorted(os.listdir(nl)):
        if f.endswith(".mps"):
            jobs.append((os.path.join(nl, f), f[:-4], [], []))
            jobs.append((os.path.join(nl, f), f[:-4], ["--no-presolve"], []))
    if "--relaxed" in sys.argv:
        d = sys.argv[sys.argv.index("--relaxed") + 1]
        for n in RELAXED:
            jobs.append((os.path.join(d, n + ".mps"), n, [], []))
    bad = bails = 0
    for path, name, extra, _ in jobs:
        ok, ho = highs_obj(path)
        if not ok:
            continue
        s, o, it, bail = ours(path, extra)
        good = s == "Optimal" and abs(o - ho) <= 1e-6 * max(1.0, abs(ho))
        bad += (not good) and name not in KNOWN
        bails += bail
        print(f"{name:22s} {'nopre' if extra else 'pre  '} {s:10s} {o:.8g} vs {ho:.8g} it={it} "
              f"{'BAIL ' if bail else ''}{'ok' if good else ('KNOWN' if name in KNOWN else 'FAIL')}", flush=True)
    print(f"FAILURES: {bad}  silent bails: {bails}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
