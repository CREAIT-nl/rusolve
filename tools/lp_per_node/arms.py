#!/usr/bin/env python3
"""Arm table: stats.work, nodes, LP iterations per arm vs the default (full-*) runs.
usage: arms.py DIR   (DIR holds full-<inst>-<n>.json and <arm>-<inst>-<n>.json)"""
import json, glob, math, os, sys
d = sys.argv[1]
rows = {}
for f in glob.glob(os.path.join(d, '*.json')):
    arm, rest = os.path.basename(f)[:-5].split('-', 1)
    inst = rest.rsplit('-', 1)[0]
    try:
        j, _ = json.JSONDecoder().raw_decode(open(f).read())
    except Exception:
        continue
    s = j.get('stats') or {}
    rows.setdefault(arm, {})[inst] = (j.get('status'), s.get('work'), s.get('nodes'), s.get('lp_iterations'))
base = rows.pop('full')
rows.pop('base', None); rows.pop('faith', None)  # node-limited decomposition runs
for arm in sorted(rows):
    print(f'\n### {arm}\n| instance | status | work / default | nodes / default | LP it / default |\n|---|---|---:|---:|---:|')
    wr = []; solved = 0; better = 0
    for inst in sorted(base):
        b = base[inst]; r = rows[arm].get(inst)
        if not r or r[1] is None:
            print(f'| {inst} | not run | | | |'); continue
        ok = r[0] == 'Optimal'; solved += ok
        w = r[1] / b[1]; wr.append(w); better += (w <= 0.9 and ok)
        print(f'| {inst} | {r[0]} | {w:.2f} | {r[2]/b[2]:.2f} | {r[3]/b[3]:.2f} |')
    if wr:
        print(f'geomean work {math.exp(sum(map(math.log, wr))/len(wr)):.3f}, optimal {solved}/{len(wr)}, instances <=0.9x work: {better}')
