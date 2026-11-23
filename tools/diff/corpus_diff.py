#!/usr/bin/env python3
"""Corpus validation gate: every manifest instance must read cleanly through
rusolve-cli AND match highspy on rows/cols/nnz; solve status+objective
agreement is checked where a solve is quick (see timeboxing below).

    python tools/diff/corpus_diff.py [--tier lpopt|mip50|all]
        [--only NAME[,NAME]] [--solve-limit 30] [--max-dense 2.5e7]
        [--on-demand] [--bin PATH] [--csv OUT]

READ (mandatory): `rusolve-cli dump <file>` exits 0 and
(n_vars, n_cons, nnz) == highspy (num_col_, num_row_, len(a_matrix_.value_)).

SOLVE (timeboxed, informational): skipped with note `skip(size)` when the
dense LP footprint m*(n+3m) exceeds --max-dense doubles (default 2.5e7 ≈ 200 MB
— our solver is dense; pds-100/irish-electricity-class models would OOM,
reading them is still mandatory). Otherwise both solvers get --solve-limit
seconds (default 30; subprocess kill for us, native time_limit for HiGHS).
Verdicts: AGREE (status match, obj rel 1e-5 when Optimal), MISMATCH (both terminal, disagree — a solver bug, reported not
crashed on), TIMEBOXED (HiGHS itself did not reach a terminal status: nothing
to compare), ERR (our side crashed/errored while HiGHS finished).

Exit 0 iff every instance READs PASS. Solver bugs surface in the table, never
block the corpus check (they are fixed separately).
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DATA = REPO / "benches" / "data"

import highspy  # noqa: E402

IORU = "InfeasibleOrUnbounded"
HIGH_TERM = {"Optimal", "Infeasible", "Unbounded", IORU}
OBJ_TOL = 1e-5


def our_dump(binary: str, path: Path, timeout: int):
    p = subprocess.run([binary, "dump", str(path)], capture_output=True,
                       text=True, timeout=timeout)
    if p.returncode != 0:
        return None, p.stderr.strip()[:200]
    out = "\n" + p.stdout  # count-based parse: supportcase10-class dumps are
    nv = _int_after(out, "\nn_vars ")   # ~20M lines — splitlines() would blow
    nc = _int_after(out, "\nn_cons ")   # up a multi-GB string list
    nnz = out.count("\nnz ")
    return (nv, nc, nnz), ""


def _int_after(text: str, marker: str):
    i = text.find(marker)
    if i < 0:
        return None
    j = text.find("\n", i + 1)
    return int(text[i + len(marker):j].split()[0])


def highs_read(path: Path):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    if h.readModel(str(path)) != highspy.HighsStatus.kOk:
        return None, None
    lp = h.getLp()
    return h, (lp.num_row_, lp.num_col_, len(lp.a_matrix_.value_))


def highs_solve(h, is_mip, limit):
    if is_mip:
        h.setOptionValue("mip_rel_gap", 0.0)
        h.setOptionValue("mip_abs_gap", 1e-9)
    h.setOptionValue("time_limit", float(limit))
    h.run()
    st = h.modelStatusToString(h.getModelStatus())
    st = {"Time limit reached": "Timeout",
          "Iteration limit reached": "IterationLimit",
          "Primal infeasible or unbounded": IORU}.get(st, st)
    obj = h.getInfo().objective_function_value if st == "Optimal" else None
    return st, obj


def our_solve(binary: str, path: Path, is_mip: bool, limit: float):
    cmd = [binary, "mip" if is_mip else "solve", str(path), "--csv"]
    if is_mip:
        cmd += ["--time-limit", str(int(limit))]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=limit + 30)
    except subprocess.TimeoutExpired:
        return "Timeout", None
    lines = p.stdout.strip().splitlines()
    if len(lines) < 2:
        return "ERR", (p.stderr.strip().splitlines() or [""])[0][:160]
    parts = lines[1].split(",", 2)
    try:
        return parts[0], float(parts[1])
    except (ValueError, IndexError):
        return "ERR", parts[0]


def solve_verdict(o_st, o_obj, h_st, h_obj):
    if h_st not in HIGH_TERM:
        return "TIMEBOXED" if o_st in ("Timeout",) else f"unref({o_st})"
    if o_st == "ERR":
        return f"ERR:{h_st}"
    st_ok = o_st == h_st or (h_st == IORU and o_st in ("Infeasible", "Unbounded"))
    if not st_ok:
        return f"MISMATCH:status({o_st} vs {h_st})"
    if h_st == "Optimal" and abs(o_obj - h_obj) > OBJ_TOL * max(1.0, abs(h_obj)):
        return f"MISMATCH:obj({o_obj:.10g} vs {h_obj:.10g})"
    return "AGREE"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tier", choices=("lpopt", "mip50", "all"), default="all")
    ap.add_argument("--only")
    ap.add_argument("--on-demand", action="store_true")
    ap.add_argument("--solve-limit", type=float, default=30.0)
    ap.add_argument("--max-dense", type=float, default=2.5e7,
                    help="skip solve above this many doubles, m*(n+3m) (default 200 MB)")
    ap.add_argument("--dump-timeout", type=int, default=600)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--csv", default=str(HERE / "corpus_diff.csv"))
    args = ap.parse_args()

    binary = args.bin or next(
        (str(p) for p in (REPO / "target/release/rusolve-cli",
                          REPO / "target/debug/rusolve-cli") if p.exists()), None)
    if binary is None:
        sys.exit("rusolve-cli not built: cargo build --release -p rusolve-cli")
    man = tomllib.loads((DATA / "manifest.toml").read_text())
    want = None if args.tier == "all" else args.tier
    only = set(args.only.split(",")) if args.only else None
    entries = [e for e in man["instance"]
               if (want is None or e["tier"] == want)
               and (only is None or e["name"] in only)
               and (args.on_demand or not e.get("on-demand"))]

    rows, all_read_ok = [], True
    for e in entries:
        t, name = e["tier"], e["name"]
        path = DATA / t / f"{name}.mps"
        if not path.exists():
            rows.append((t, name, "-", "-", "-", "FAIL", "missing",
                         "not downloaded: see benches/data/ATTRIBUTION.md"))
            all_read_ok = False
            continue
        try:
            ours, err = our_dump(binary, path, args.dump_timeout)
        except subprocess.TimeoutExpired:
            ours, err = None, "dump timeout"
        h, hcounts = highs_read(path)
        if h is None:
            rows.append((t, name, "-", "-", "-", "FAIL", "ORACLE-READ-FAIL",
                         "highspy readModel failed"))
            all_read_ok = False
            continue
        if ours is None or ours != (hcounts[1], hcounts[0], hcounts[2]):
            rows.append((t, name,
                         "-" if ours is None else f"{ours[1]}/{hcounts[0]}",
                         "-" if ours is None else f"{ours[0]}/{hcounts[1]}",
                         "-" if ours is None else f"{ours[2]}/{hcounts[2]}",
                         "FAIL", "READ", err))
            all_read_ok = False
            continue
        # ---- solve stage (informational)
        m, n, _ = hcounts
        if m * (n + 3 * m) > args.max_dense:
            rows.append((t, name, f"{ours[1]}/{m}", f"{ours[0]}/{n}",
                         f"{ours[2]}/{_}", "PASS", "skip(size)",
                         f"dense {m*(n+3*m):.1e} > {args.max_dense:.0e}"))
            continue
        lp = h.getLp()
        is_mip = any(v != highspy.HighsVarType.kContinuous for v in lp.integrality_)
        h_st, h_obj = highs_solve(h, is_mip, args.solve_limit)
        o_st, o_obj = our_solve(binary, path, is_mip, args.solve_limit)
        note = o_obj if isinstance(o_obj, str) else ""
        v = solve_verdict(o_st, o_obj if not isinstance(o_obj, str) else None,
                          h_st, h_obj)
        rows.append((t, name, f"{ours[1]}/{m}", f"{ours[0]}/{n}", f"{ours[2]}/{_}",
                     "PASS", v, note))

    print(f"{'tier':6s} {'instance':26s} {'rows u/h':>16s} {'cols u/h':>16s} "
          f"{'nnz u/h':>18s} {'READ':5s} {'SOLVE':22s}")
    for r in rows:
        print(f"{r[0]:6s} {r[1]:26s} {r[2]:>16s} {r[3]:>16s} {r[4]:>18s} "
              f"{r[5]:5s} {r[6]:22s}" + (f"  {r[7][:60]}" if r[7] else ""))
    with open(args.csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tier", "instance", "rows", "cols", "nnz", "read", "solve", "note"])
        w.writerows(rows)
    n_pass = sum(1 for r in rows if r[5] == "PASS")
    print(f"\nread: {n_pass}/{len(rows)} PASS   solve: "
          + "  ".join(f"{v}={sum(1 for r in rows if r[6].split('(')[0].split(':')[0] == v)}"
                      for v in ("AGREE", "MISMATCH", "TIMEBOXED", "ERR", "skip")))
    print(f"table -> {args.csv}")
    print("ALL READ PASS" if all_read_ok else "READ FAILURES PRESENT")
    return 0 if all_read_ok else 1


if __name__ == "__main__":
    sys.exit(main())
