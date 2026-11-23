#!/usr/bin/env python3
"""IPM benchmark: the rusolve interior-point solver against its dual simplex and HiGHS IPM.

For each instance runs, one process at a time, with a peak-RSS watchdog:
  * rusolve IPM        `rusolve-cli solve <f> --ipm`
  * rusolve dual simplex `rusolve-cli solve <f>`        (unless --no-simplex)
  * HiGHS ipm           highspy, solver=ipm, threads=1, run_crossover=off
                        (MIP files: solve_relaxation=true)  (unless --no-highs)
Records status, objective, iterations, wall, peak RSS and (rusolve IPM) the
Cholesky statistics from the `rusolve ipm stats:` stderr line. Output: JSON
to --out (never into the repo checkout's tracked files).

usage: python tools/bench/ipm_bench.py --bin target/release/rusolve-cli --out OUT.json
                                    [--timeout S] [--rss-kill-gb 16] files...
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

HIGHS_SNIPPET = r"""
import json, sys, time, highspy
f, tl = sys.argv[1], float(sys.argv[2])
h = highspy.Highs()
h.setOptionValue("output_flag", False)
h.setOptionValue("threads", 1)
h.setOptionValue("solver", "ipm")
h.setOptionValue("run_crossover", "off")
h.setOptionValue("time_limit", tl)
h.setOptionValue("solve_relaxation", True)
t0 = time.perf_counter()
h.readModel(f)
h.run()
t = time.perf_counter() - t0
info = h.getInfo()
print(json.dumps({"status": h.modelStatusToString(h.getModelStatus()),
                  "obj": info.objective_function_value,
                  "iters": info.ipm_iteration_count, "time": t}))
"""


def tree_rss_kb(pid):
    """RSS (KiB) of pid's children (the /usr/bin/time payload)."""
    kids = subprocess.run(["pgrep", "-P", str(pid)], capture_output=True, text=True).stdout.split()
    tot = 0
    for k in kids:
        try:
            tot += int(subprocess.run(["ps", "-o", "rss=", "-p", k],
                                      capture_output=True, text=True).stdout.strip() or 0)
        except ValueError:
            pass
    return tot


def run_watched(cmd, timeout, rss_kill_kb):
    """Run cmd under /usr/bin/time -l (exact peak RSS), poll RSS every
    0.2 s as a watchdog, kill above rss_kill_kb."""
    t0 = time.perf_counter()
    p = subprocess.Popen(["/usr/bin/time", "-l"] + cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True)
    peak = 0
    killed = None
    while True:
        try:
            out, err = p.communicate(timeout=0.2)
            break
        except subprocess.TimeoutExpired:
            pass
        rss = tree_rss_kb(p.pid)
        peak = max(peak, rss)
        if rss > rss_kill_kb or time.perf_counter() - t0 > timeout + 30:
            killed = "rss" if rss > rss_kill_kb else "timeout"
            subprocess.run(["pkill", "-9", "-P", str(p.pid)])
            p.kill()
    wall = time.perf_counter() - t0
    m = re.search(r"(\d+)\s+maximum resident set size", err or "")
    if m:
        peak = max(peak, int(m.group(1)) // 1024)
    return out, err, wall, peak, killed, p.returncode


def parse_rusolve(out, err):
    r = {}
    m = re.search(r"^status:\s+(\S+)", out, re.M)
    r["status"] = m.group(1) if m else "ERROR"
    m = re.search(r"^objective: (\S+)", out, re.M)
    r["obj"] = float(m.group(1)) if m else None
    m = re.search(r"^iterations: (\d+) \(phase1\) \+ (\d+) \(phase2\)", out, re.M)
    r["iters"] = int(m.group(1)) + int(m.group(2)) if m else None
    m = re.search(r"rusolve ipm stats: (.*)", err)
    if m:
        s = m.group(1)
        for key, pat in [("m", r"\bm (\d+)"), ("n", r"\bn (\d+)"), ("nnz_a", r"nnz\(A\) (\d+)"),
                         ("nnz_m", r"nnz\(M\) (\d+)"), ("nnz_l", r"nnz\(L\) (\d+)"),
                         ("supernodes", r"supernodes (\d+)"), ("dense_cols", r"dense_cols (\d+)"), ("dense_split", r"dense_split (\d+)"),
                         ("correctors", r"correctors (\d+)"),
                         ("skipped", r"skipped_pivots (\d+)")]:
            mm = re.search(pat, s)
            if mm:
                r[key] = int(mm.group(1))
        for key in ["order", "symbolic", "factor", "solve", "ipm"]:
            mm = re.search(key + r" ([\d.]+)s", s)
            if mm:
                r["t_" + key] = float(mm.group(1))
        r["retried"] = "retried true" in s
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--rss-kill-gb", type=float, default=16.0)
    ap.add_argument("--no-simplex", action="store_true")
    ap.add_argument("--no-highs", action="store_true")
    ap.add_argument("--ipm-flags", default="")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    kill_kb = int(a.rss_kill_gb * 1024 * 1024)
    res = {}
    if os.path.exists(a.out):
        res = json.load(open(a.out))
    for f in a.files:
        name = os.path.basename(f).replace(".mps", "").replace(".gz", "")
        row = res.get(name, {})
        tl = str(a.timeout)
        out, err, wall, peak, killed, rc = run_watched(
            [a.bin, "solve", f, "--ipm", "--time-limit", tl] + a.ipm_flags.split(),
            a.timeout, kill_kb)
        ipm = parse_rusolve(out, err)
        ipm.update(wall=wall, rss_mb=peak / 1024, killed=killed)
        row["ipm"] = ipm
        if not a.no_simplex:
            out, err, wall, peak, killed, rc = run_watched(
                [a.bin, "solve", f, "--time-limit", tl], a.timeout, kill_kb)
            sx = parse_rusolve(out, err)
            sx.update(wall=wall, rss_mb=peak / 1024, killed=killed)
            row["simplex"] = sx
        if not a.no_highs:
            out, err, wall, peak, killed, rc = run_watched(
                [sys.executable, "-c", HIGHS_SNIPPET, f, tl], a.timeout, kill_kb)
            try:
                hs = json.loads(out.strip().splitlines()[-1])
            except (ValueError, IndexError):
                hs = {"status": "ERROR", "err": err[-300:]}
            hs.update(wall=wall, rss_mb=peak / 1024, killed=killed)
            row["highs_ipm"] = hs
        res[name] = row
        json.dump(res, open(a.out, "w"), indent=1, sort_keys=True)
        print(name, json.dumps({k: (v.get("status"), v.get("obj"), v.get("iters"),
                                    round(v.get("wall", 0), 3), round(v.get("rss_mb", 0)))
                                for k, v in row.items()}), flush=True)


if __name__ == "__main__":
    main()
