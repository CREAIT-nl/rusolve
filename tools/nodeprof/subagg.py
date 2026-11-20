#!/usr/bin/env python3
"""Mean subsystem share (percent of run) port vs C++ over instances; ratio of shares.
usage: subagg.py inst:nodes ...   (needs $NODEPROF_DIR/<arm>-<inst>-<nodes>.json from cmp.py, default out/)"""
import json, sys, collections, math, os
specs = sys.argv[1:]
sh = {arm: collections.defaultdict(list) for arm in ('port', 'cpp')}
for arm in sh:
    for sp in specs:
        i, n = sp.split(':')
        d = json.load(open(f"{os.environ.get('NODEPROF_DIR', 'out')}/{arm}-{i}-{n}.json"))
        for k, v in d['sub'].items(): sh[arm][k].append(100 * v / d['n'])
        for k in set(sh[arm]) - set(d['sub']): sh[arm][k].append(0)
print(f'{"subsystem":18s} {"port%":>7s} {"cpp%":>7s} {"ratio":>6s}   per-instance port/cpp ratio')
for k in sorted(set(sh['port']) | set(sh['cpp'])):
    p = sh['port'][k] + [0] * (len(specs) - len(sh['port'][k])); c = sh['cpp'][k] + [0] * (len(specs) - len(sh['cpp'][k]))
    mp, mc = sum(p) / len(p), sum(c) / len(c)
    if mp + mc < 0.3: continue
    per = ' '.join(f'{(a / b if b > 0.05 else float("nan")):5.2f}' for a, b in zip(p, c))
    print(f'{k:18s} {mp:7.2f} {mc:7.2f} {mp / mc if mc else float("inf"):6.2f}   {per}')
