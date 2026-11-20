#!/usr/bin/env python3
"""Runs the fenced ```python blocks of docs/skills/*/SKILL.md.

A playbook's code is only worth shipping if it runs and gives the answer the text states, so
every block is an executable program that ends in `assert`s on the optimum. This script extracts
the blocks and runs each in a fresh interpreter, with the rusolve Python package importable:

    python3 evals/skill_code.py                     # every skill
    python3 evals/skill_code.py lot-sizing          # one skill (directory name under docs/skills)

Needs `import rusolve` (maturin develop -m crates/rusolve-py/Cargo.toml). Exit 1 when a block
fails. A block fenced as ```python no-run is shown to the reader and skipped here. Standard
library only.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKILLS = os.path.join(ROOT, "docs", "skills")
FENCE = re.compile(r"^```python([^\n]*)\n(.*?)^```", re.S | re.M)


def blocks(path):
    text = open(path).read()
    return [(m.group(1).strip(), m.group(2), text[:m.start()].count("\n") + 2) for m in FENCE.finditer(text)]


def main():
    want = set(a for a in sys.argv[1:] if not a.startswith("-"))
    failed = ran = 0
    for d in sorted(os.listdir(SKILLS)):
        md = os.path.join(SKILLS, d, "SKILL.md")
        if not os.path.isfile(md) or (want and d not in want):
            continue
        for flag, code, line in blocks(md):
            if flag == "no-run":
                continue
            ran += 1
            r = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, timeout=300)
            ok = r.returncode == 0
            failed += not ok
            print(f"{'ok  ' if ok else 'FAIL'} {d}/SKILL.md:{line}")
            if r.stdout.strip():
                print("     " + r.stdout.strip().replace("\n", "\n     "))
            if not ok:
                print("     " + r.stderr.strip()[-800:].replace("\n", "\n     "))
    print(f"{ran - failed}/{ran} blocks ok")
    return 1 if failed or not ran else 0


if __name__ == "__main__":
    sys.exit(main())
