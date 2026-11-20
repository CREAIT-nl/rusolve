#!/usr/bin/env python3
"""Generate the multi-knapsack instance family that `rusolve tune` is tested on, a
line-for-line port of `knapsack()` in crates/rusolve-cli/tests/tune.rs (an LCG with Numerical Recipes constants, Press et al., "Numerical Recipes" 2nd ed. 7.1).
  python3 tools/gen_tune_family.py DIR [--count 24] [--n 40] [--m 5]
writes DIR/k00.mps .. k{count-1}.mps (seed = file index)."""
import argparse
import os


def knapsack(seed, n, m):
    x = 12345 + seed * 7919

    def r():
        nonlocal x
        x = (x * 1103515245 + 12345) % (1 << 31)
        return x >> 8

    w = [[r() % 47 + 10 for _ in range(n)] for _ in range(m)]
    s = "NAME KNAP\nROWS\n N  cost\n"
    for i in range(m):
        s += f" L  k{i}\n"
    s += "COLUMNS\n    MARKER  'MARKER'  'INTORG'\n"
    for j in range(n):
        s += f"    x{j}  cost  {-(r() % 40 + 10)}\n"
        for i, row in enumerate(w):
            s += f"    x{j}  k{i}  {row[j]}\n"
    s += "    MARKER  'MARKER'  'INTEND'\nRHS\n"
    for i, row in enumerate(w):
        s += f"    rhs  k{i}  {sum(row) // 3}\n"
    s += "BOUNDS\n"
    for j in range(n):
        s += f" UP bnd x{j} 1\n"
    return s + "ENDATA\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--count", type=int, default=24)
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--m", type=int, default=5)
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    for seed in range(a.count):
        with open(os.path.join(a.dir, f"k{seed:02d}.mps"), "w") as f:
            f.write(knapsack(seed, a.n, a.m))
    print(f"wrote {a.count} instances to {a.dir}")


if __name__ == "__main__":
    main()
