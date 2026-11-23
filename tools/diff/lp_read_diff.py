#!/usr/bin/env python3
"""LP-reader acceptance gate: rusolve-io LP parser vs highspy readModel.

Corpus: tools/diff/data/lp_samples/*.lp — handcrafted small LPs covering the
shared CPLEX-LP-format grammar (objective label/constant, named and unnamed
rows, L/G/E senses, implicit coefficients, scientific notation, Bounds with
free / one-sided / equality / range chains, no-Bounds defaults, comments,
punctuated names, wrapped statements, case-insensitive keywords, Binary and
General sections, unbounded and infeasible models).

Deliberately OUT of the corpus (our reader supports them, HiGHS's LP parser
does not — unit-tested in crates/rusolve-io/src/lp.rs instead): two-sided
constraint rows (`lo <= expr <= hi`), `*` coefficients, `in [l, u]` bounds,
comma-separated bounds lists, `Free` sections, GLPK `default` bounds,
space-containing names, bare-expression rows.

For every sample we check:
  (1) row/col counts agree with highspy readModel (`dump` vs lp.num_row_/
      num_col_);
  (2) after solving, status agrees and, when Optimal, objectives agree to
      rel 1e-6 (`solve --csv` for LPs, `mip --csv` for integer models);
  (3) transparent gzip: our dump of <sample>.lp.gz and of a bare .gz copy
      (content-sniffed) equals the dump of the plain file.

Usage: python tools/diff/lp_read_diff.py
Exit 0 iff ALL PASS.
"""
import gzip
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SAMPLES = os.path.join(HERE, "data", "lp_samples")
BIN = os.path.join(REPO, "target", "debug", "rusolve-cli")

import highspy
from highspy import HighsModelStatus, HighsStatus, HighsVarType


def dump(path):
    """`rusolve-cli dump <path>` -> (n_vars, n_cons, raw lines)."""
    p = subprocess.run([BIN, "dump", path], capture_output=True, text=True)
    if p.returncode != 0:
        return None, None, p.stderr
    nv = nc = None
    for line in p.stdout.splitlines():
        if line.startswith("n_vars "):
            nv = int(line.split()[1])
        elif line.startswith("n_cons "):
            nc = int(line.split()[1])
    return nv, nc, p.stdout


def solve_ours(path, is_mip):
    if is_mip:
        cmd = [BIN, "mip", path, "--csv", "--node-limit", "20000",
               "--time-limit", "20"]
    else:
        cmd = [BIN, "solve", path, "--csv"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    lines = p.stdout.strip().splitlines()
    if len(lines) < 2:
        return "ERR", None, p.stderr
    parts = lines[1].split(",")
    return parts[0], float(parts[1]), p.stderr


def solve_highs(h):
    h.run()
    st = h.getModelStatus()
    obj = h.getObjectiveValue() if st == HighsModelStatus.kOptimal else None
    # HiGHS "kOptimal" with |obj| >= 1e15 is unbounded-direction runaway
    # (same guard as simplex_diff.py).
    if st == HighsModelStatus.kOptimal and abs(obj) >= 1e15:
        return HighsModelStatus.kUnbounded, None
    return st, obj


def verdict(our_status, our_obj, h_status, h_obj):
    if h_status == HighsModelStatus.kOptimal:
        return (our_status == "Optimal"
                and abs(our_obj - h_obj) <= 1e-6 * max(1.0, abs(h_obj)))
    if h_status == HighsModelStatus.kInfeasible:
        return our_status == "Infeasible"
    if h_status == HighsModelStatus.kUnbounded:
        return our_status == "Unbounded"
    if h_status == HighsModelStatus.kUnboundedOrInfeasible:
        return our_status in ("Infeasible", "Unbounded")
    return False


def main():
    files = sorted(f for f in os.listdir(SAMPLES) if f.endswith(".lp"))
    if not files:
        print(f"no samples under {SAMPLES}")
        sys.exit(1)
    rows = []
    all_pass = True
    with tempfile.TemporaryDirectory(prefix="lp_read_diff.") as tmp:
        for fname in files:
            path = os.path.join(SAMPLES, fname)
            name = fname[:-3]

            h = highspy.Highs()
            h.setOptionValue("output_flag", False)
            h.setOptionValue("threads", 1)
            h.setOptionValue("mip_rel_gap", 0.0)
            h.setOptionValue("time_limit", 30.0)
            read = h.readModel(path)
            if read != HighsStatus.kOk:
                rows.append((name, "-", "-", "ORACLE-READ-FAIL", "-", "FAIL"))
                all_pass = False
                continue
            lp = h.getLp()
            is_mip = any(v != HighsVarType.kContinuous
                         for v in lp.integrality_)

            nv, nc, dump_raw = dump(path)
            counts_ok = (nv, nc) == (lp.num_col_, lp.num_row_)

            h_status, h_obj = solve_highs(h)
            o_status, o_obj, err = solve_ours(path, is_mip)
            solve_ok = counts_ok and verdict(o_status, o_obj, h_status, h_obj)

            # transparent .gz: plain .lp.gz and bare .gz (sniffed) must
            # dump identically to the plain file.
            gz_path = os.path.join(tmp, name + ".lp.gz")
            sniff_path = os.path.join(tmp, name + ".gz")
            with open(path, "rb") as f, gzip.open(gz_path, "wb") as g, \
                    gzip.open(sniff_path, "wb") as s:
                data = f.read()
                g.write(data)
                s.write(data)
            gz_ok = (dump(gz_path)[2] == dump_raw) and \
                    (dump(sniff_path)[2] == dump_raw)

            ok = solve_ok and gz_ok
            all_pass &= ok
            rows.append((
                name,
                f"{nc}/{lp.num_row_}" if nv is not None else "-",
                f"{nv}/{lp.num_col_}" if nv is not None else "-",
                o_status,
                f"{'%.6g' % o_obj if o_obj is not None else '-'}/"
                f"{'%.6g' % h_obj if h_obj is not None else str(h_status).replace('HighsModelStatus.', '')}",
                ("PASS" if ok else f"FAIL{' counts' if not counts_ok else ''}"
                                   f"{' gz' if not gz_ok else ''}"),
            ))

    print(f"{'sample':>22} {'rows':>8} {'cols':>8} {'status':>12} "
          f"{'obj ours/highs':>28}")
    for r in rows:
        print(f"{r[0]:>22} {r[1]:>8} {r[2]:>8} {r[3]:>12} {r[4]:>28}  {r[5]}")
    print("ALL PASS" if all_pass else "FAILURES PRESENT")
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
