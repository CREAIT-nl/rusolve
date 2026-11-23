"""Write infeasible / maximisation variants of netlib MPS files for the
simplex traces, so the fixtures reach the non-optimal exits that the stock
netlib instances (all Optimal) never do.

usage: python3 -I mkvariants.py <netlib_dir> <out_dir>

For each LP in ZSTAR:
  <lp>max.mps  OBJSENSE MAX inserted after NAME. Some are unbounded, and
               the dual run then hands over to primal.
  <lp>cut.mps  an extra row ZCUT: c'x <= z* - (0.01|z*| + 1), with c the
               objective row. It cuts off every feasible point, so the LP is
               primal infeasible (dual ray proof / primal phase-1 exit).

z* is the optimal value published with the netlib LP set
(netlib.org/lp/data/readme). The output is a pure function of the inputs.
"""
import os
import sys

ZSTAR = {
    "afiro": -464.7531428571,
    "adlittle": 225494.96316,
    "sc50a": -64.575077059,
    "sc105": -52.202061212,
    "blend": -30.812149846,
    "kb2": -1749.9001299,
    "share2b": -415.73224074,
}


def max_variant(lines):
    return [lines[0], "OBJSENSE", "    MAX"] + lines[1:]


def cut_variant(lines, z):
    section = None
    obj = None
    rhs_set = None
    rhs_written = False
    res = []
    for ln in lines:
        if ln and not ln[0].isspace():
            head = ln.split()[0]
            if head == "COLUMNS":
                res.append(" L  ZCUT")  # close ROWS with the cut row
            if section == "RHS" and head != "RHS":
                res.append("    %-8s  %-8s  %.12g" % (rhs_set or "RHS", "ZCUT", z - (0.01 * abs(z) + 1.0)))
                rhs_written = True
            section = head
            res.append(ln)
            continue
        tok = ln.split()
        if section == "ROWS" and tok and tok[0] == "N" and obj is None:
            obj = tok[1]
        if section == "COLUMNS" and len(tok) >= 3 and "MARKER" not in ln:
            res.append(ln)
            for r, v in zip(tok[1::2], tok[2::2]):
                if r == obj:
                    res.append("    %-8s  %-8s  %s" % (tok[0], "ZCUT", v))
            continue
        if section == "RHS" and tok and rhs_set is None and len(tok) % 2 == 1:
            rhs_set = tok[0]
        res.append(ln)
    if obj is None or not rhs_written:
        sys.exit("mkvariants: could not place the ZCUT row/rhs")
    return res


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    for lp, z in ZSTAR.items():
        with open(os.path.join(src, lp + ".mps")) as f:
            lines = f.read().splitlines()
        for suffix, body in (("max", max_variant(lines)), ("cut", cut_variant(lines, z))):
            with open(os.path.join(out, lp + suffix + ".mps"), "w") as f:
                f.write("\n".join(body) + "\n")


if __name__ == "__main__":
    main()
