#!/usr/bin/env python3
"""Independent MIP solution checker + MPS writer for the mip_diff gate.

Deliberately self-contained: the checker re-reads the MPS FILE with its own
small free-format reader (no rusolve code, no solver library), so a verdict
never depends on how our solver or an oracle interpreted the model.

Semantics (MPS standard as implemented by rusolve-io / HiGHS / SCIP):
  * default column bounds [0, +inf); columns inside MARKER INTORG/INTEND are
    integer and, when no BOUNDS line names them, BINARY [0, 1] (HiGHS
    kintegerVarsInColumnsAreBinary convention, adopted by rusolve-io in
    5ee572f); any BOUNDS line clears that default. BV = binary, LI/UI =
    integer bound;
  * RANGES: L -> [rhs-|R|, rhs], G -> [rhs, rhs+|R|],
            E -> [rhs, rhs+R] if R > 0 else [rhs+R, rhs];
  * an RHS entry on the objective row is a constant offset of -value.

Tolerances (MIP_diff contract): integrality |x - round(x)| <= 1e-6; bound
and row feasibility 1e-6 scaled by max(1, |bound|, largest |a_ij x_j| term)
so big-M rows are judged relative to their own magnitude; objective
recomputed from x and compared to the reported value at 1e-6 relative to
max(1, sum |c_j x_j|).

Usage: python tools/diff/mip_check.py <model.mps> <rusolve-cli mip --csv output file>
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

INF = float("inf")
INT_TOL = 1e-6
FEAS_TOL = 1e-6
OBJ_TOL = 1e-6


# --------------------------------------------------------------------------- model


@dataclass
class MpsModel:
    """Column/row data exactly as the MPS file states it."""

    name: str = ""
    maximize: bool = False
    colnames: List[str] = field(default_factory=list)
    lb: List[float] = field(default_factory=list)
    ub: List[float] = field(default_factory=list)
    is_int: List[bool] = field(default_factory=list)
    obj: List[float] = field(default_factory=list)
    obj_offset: float = 0.0
    rownames: List[str] = field(default_factory=list)
    rows: List[List[Tuple[int, float]]] = field(default_factory=list)
    row_lo: List[float] = field(default_factory=list)
    row_hi: List[float] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.colnames)

    @property
    def m(self) -> int:
        return len(self.rownames)


def row_bounds(sense: str, rhs: float, rng: float) -> Tuple[float, float]:
    """(lo, hi) of one MPS row, RANGES semantics as in simplex_diff."""
    if rng == 0.0:
        return {"L": (-INF, rhs), "G": (rhs, INF), "E": (rhs, rhs)}[sense]
    if sense == "L":
        return (rhs - abs(rng), rhs)
    if sense == "G":
        return (rhs, rhs + abs(rng))
    return (rhs, rhs + rng) if rng > 0 else (rhs + rng, rhs)


def read_mps(path: str) -> MpsModel:
    """Free-format MPS reader (fixed-format files without spaces in names too)."""
    mdl = MpsModel()
    section = None
    obj_row: Optional[str] = None
    row_idx: Dict[str, int] = {}
    senses: List[str] = []
    col_idx: Dict[str, int] = {}
    rhs: Dict[int, float] = {}
    ranges: Dict[int, float] = {}
    in_int = False
    has_bound: Dict[int, bool] = {}
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("*"):
                continue
            t = line.split()
            if not line[0].isspace():
                head = t[0].upper()
                if head == "NAME":
                    mdl.name = t[1] if len(t) > 1 else ""
                    continue
                if head == "OBJSENSE":
                    section = "OBJSENSE"
                    if len(t) > 1:
                        mdl.maximize = t[1].upper() in ("MAX", "MAXIMIZE")
                        section = None
                    continue
                if head in ("ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS"):
                    section = head
                    continue
                if head == "ENDATA":
                    break
                raise ValueError(f"{path}: unknown section {t[0]!r}")
            if section == "OBJSENSE":
                mdl.maximize = t[0].upper() in ("MAX", "MAXIMIZE")
            elif section == "ROWS":
                s, name = t[0].upper(), t[1]
                if s == "N":
                    if obj_row is None:
                        obj_row = name
                    continue
                row_idx[name] = len(mdl.rownames)
                mdl.rownames.append(name)
                senses.append(s)
                mdl.rows.append([])
            elif section == "COLUMNS":
                if len(t) >= 3 and t[1].strip("'").upper() == "MARKER":
                    in_int = t[2].strip("'").upper() == "INTORG"
                    continue
                c = t[0]
                if c not in col_idx:
                    col_idx[c] = len(mdl.colnames)
                    mdl.colnames.append(c)
                    mdl.lb.append(0.0)
                    mdl.ub.append(INF)
                    mdl.is_int.append(in_int)
                    mdl.obj.append(0.0)
                j = col_idx[c]
                for k in range(1, len(t) - 1, 2):
                    r, v = t[k], float(t[k + 1])
                    if r == obj_row:
                        mdl.obj[j] += v
                    elif r in row_idx:
                        mdl.rows[row_idx[r]].append((j, v))
                    else:
                        raise ValueError(f"{path}: column {c} references unknown row {r}")
            elif section in ("RHS", "RANGES"):
                pairs = t[1:] if len(t) % 2 == 1 else t
                for k in range(0, len(pairs) - 1, 2):
                    r, v = pairs[k], float(pairs[k + 1])
                    if section == "RHS" and r == obj_row:
                        mdl.obj_offset = -v
                    elif r in row_idx:
                        (rhs if section == "RHS" else ranges)[row_idx[r]] = v
                    else:
                        raise ValueError(f"{path}: {section} references unknown row {r}")
            elif section == "BOUNDS":
                ty = t[0].upper()
                # TYPE [setname] COL [value]
                if t[1] in col_idx and (len(t) == 2 or t[2] not in col_idx):
                    c, val = t[1], (float(t[2]) if len(t) > 2 else None)
                else:
                    c, val = t[2], (float(t[3]) if len(t) > 3 else None)
                j = col_idx[c]
                has_bound[j] = True
                if ty == "UP":
                    mdl.ub[j] = val
                elif ty == "LO":
                    mdl.lb[j] = val
                elif ty == "FX":
                    mdl.lb[j] = mdl.ub[j] = val
                elif ty == "FR":
                    mdl.lb[j], mdl.ub[j] = -INF, INF
                elif ty == "MI":
                    mdl.lb[j] = -INF
                elif ty == "PL":
                    mdl.ub[j] = INF
                elif ty == "BV":
                    mdl.lb[j], mdl.ub[j], mdl.is_int[j] = 0.0, 1.0, True
                elif ty == "LI":
                    mdl.lb[j], mdl.is_int[j] = val, True
                elif ty == "UI":
                    mdl.ub[j], mdl.is_int[j] = val, True
                else:
                    raise ValueError(f"{path}: unsupported bound type {ty}")
    for j in range(mdl.n):
        if mdl.is_int[j] and not has_bound.get(j):
            mdl.ub[j] = 1.0
    for i, s in enumerate(senses):
        lo, hi = row_bounds(s, rhs.get(i, 0.0), ranges.get(i, 0.0))
        mdl.row_lo.append(lo)
        mdl.row_hi.append(hi)
    return mdl


# --------------------------------------------------------------------------- writer


def write_mps(case: dict, path: str) -> None:
    """Write a case dict (see mip_diff.py) as free MPS.

    Integer columns are wrapped in MARKER INTORG/INTEND blocks; binaries are
    written as BV when case['bv'] is true, else as marker + UP 1. Every
    integer column gets at least one BOUNDS line (LO 0 when nothing else is
    needed), so the marker-means-binary default never applies implicitly.
    """
    n, m = case["n"], case["m"]
    col_rows: List[List[Tuple[int, float]]] = [[] for _ in range(n)]
    for i, row in enumerate(case["rows"]):
        for j, a in row:
            if a != 0.0:
                col_rows[j].append((i, a))
    L = [f"NAME {case.get('name', 'CASE')}",
         f"OBJSENSE\n    {'MAX' if case['maximize'] else 'MIN'}",
         "ROWS", " N  OBJ"]
    L += [f" {case['senses'][i]}  R{i}" for i in range(m)]
    L.append("COLUMNS")
    in_int = False
    for j in range(n):
        if case["is_int"][j] != in_int:
            L.append(f"    M{j}  'MARKER'  '{'INTORG' if case['is_int'][j] else 'INTEND'}'")
            in_int = case["is_int"][j]
        entries = []
        if case["obj"][j] != 0.0 or not col_rows[j]:
            entries.append(f"    X{j}  OBJ  {case['obj'][j]!r}")
        entries += [f"    X{j}  R{i}  {a!r}" for i, a in col_rows[j]]
        L += entries
    if in_int:
        L.append(f"    M{n}  'MARKER'  'INTEND'")
    L.append("RHS")
    L += [f"    RHS  R{i}  {case['rhs'][i]!r}" for i in range(m) if case["rhs"][i] != 0.0]
    if any(case["ranges"]):
        L.append("RANGES")
        L += [f"    RNG  R{i}  {case['ranges'][i]!r}" for i in range(m) if case["ranges"][i]]
    L.append("BOUNDS")
    for j in range(n):
        lo, hi = case["lb"][j], case["ub"][j]
        if case["is_int"][j] and case.get("bv") and lo == 0.0 and hi == 1.0:
            L.append(f" BV BND X{j}")
        elif lo == -INF and hi == INF:
            L.append(f" FR BND X{j}")
        elif lo == hi:
            L.append(f" FX BND X{j} {lo!r}")
        else:
            if lo == -INF:
                L.append(f" MI BND X{j}")
            elif lo != 0.0:
                L.append(f" LO BND X{j} {lo!r}")
            if hi != INF:
                L.append(f" UP BND X{j} {hi!r}")
            if case["is_int"][j] and lo == 0.0 and hi == INF:
                L.append(f" LO BND X{j} 0.0")
    L.append("ENDATA")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


# --------------------------------------------------------------------------- checker


@dataclass
class CheckResult:
    ok: bool
    obj: float  # objective recomputed from x (incl. offset)
    problems: List[str]
    max_int_viol: float = 0.0
    max_bound_viol: float = 0.0
    max_row_viol: float = 0.0  # scaled

    def short(self) -> str:
        if self.ok:
            return "ok"
        return self.problems[0] + (f" (+{len(self.problems) - 1})" if len(self.problems) > 1 else "")


def check_solution(mdl: MpsModel, x: List[float], reported_obj: Optional[float] = None,
                   feas_tol: float = FEAS_TOL, int_tol: float = INT_TOL,
                   obj_tol: float = OBJ_TOL) -> CheckResult:
    """Feasibility + integrality + objective of x against the MPS-level model."""
    probs: List[str] = []
    if len(x) != mdl.n:
        return CheckResult(False, math.nan, [f"len(x)={len(x)} != n={mdl.n}"])
    if any(not math.isfinite(v) for v in x):
        return CheckResult(False, math.nan, ["non-finite x entry"])
    res = CheckResult(True, 0.0, probs)
    for j, v in enumerate(x):
        lo, hi = mdl.lb[j], mdl.ub[j]
        viol = max(lo - v, v - hi, 0.0)
        scale = max(1.0, abs(lo) if lo > -INF else 0.0, abs(hi) if hi < INF else 0.0)
        res.max_bound_viol = max(res.max_bound_viol, viol / scale)
        if viol > feas_tol * scale:
            probs.append(f"bound {mdl.colnames[j]}={v!r} not in [{lo}, {hi}]")
        if mdl.is_int[j]:
            iv = abs(v - round(v))
            res.max_int_viol = max(res.max_int_viol, iv)
            if iv > int_tol:
                probs.append(f"integrality {mdl.colnames[j]}={v!r}")
    for i, row in enumerate(mdl.rows):
        act = math.fsum(a * x[j] for j, a in row)
        lo, hi = mdl.row_lo[i], mdl.row_hi[i]
        biggest = max((abs(a * x[j]) for j, a in row), default=0.0)
        scale = max(1.0, biggest, abs(lo) if lo > -INF else 0.0, abs(hi) if hi < INF else 0.0)
        viol = max(lo - act, act - hi, 0.0)
        res.max_row_viol = max(res.max_row_viol, viol / scale)
        if viol > feas_tol * scale:
            probs.append(f"row {mdl.rownames[i]} act={act!r} not in [{lo}, {hi}]")
    terms = [c * v for c, v in zip(mdl.obj, x)]
    res.obj = math.fsum(terms) + mdl.obj_offset
    if reported_obj is not None:
        scale = max(1.0, math.fsum(abs(t) for t in terms), abs(reported_obj))
        if not math.isfinite(reported_obj) or abs(res.obj - reported_obj) > obj_tol * scale:
            probs.append(f"objective: recomputed {res.obj!r} != reported {reported_obj!r}")
    res.ok = not probs
    return res


def parse_vector(s: str) -> List[float]:
    s = s.strip().strip("[]").strip()
    return [float(t) for t in s.split(";")] if s else []


def main(argv: List[str]) -> int:
    if len(argv) != 3:
        print(__doc__)
        return 2
    mdl = read_mps(argv[1])
    with open(argv[2]) as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    if lines and lines[0].startswith("status"):
        lines = lines[1:]
    parts = lines[0].split(",", 6)
    status, obj, x = parts[0], float(parts[1]), parse_vector(parts[6])
    r = check_solution(mdl, x, obj)
    print(f"status={status} reported_obj={obj!r} recomputed_obj={r.obj!r}")
    print(f"max viol: int {r.max_int_viol:.2e}  bound {r.max_bound_viol:.2e}  "
          f"row(scaled) {r.max_row_viol:.2e}")
    for p in r.problems[:20]:
        print("  " + p)
    print("FEASIBLE" if r.ok else "INFEASIBLE / MISMATCH")
    return 0 if r.ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
