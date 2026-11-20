#!/usr/bin/env python3
"""Function-level port-vs-C++ share gaps, averaged over instances.
usage: gaps.py inst:nodes ... [--min 0.15] [--key incl|self]
Each run's sample counts become shares of that run's total, so uniform machine-load noise
(both arms scale together) drops out; the table is the mean share over the instances.
Names are normalised to class.method (lower case, no underscores) so Rust `<..::HFactor>::ftran_u`
and C++ `HFactor::ftranU(...)` meet. Run from the checkout root after cmp.py
(cmp.py writes $NODEPROF_DIR/<arm>-<inst>-<nodes>.json, default out/)."""
import json, re, sys, collections, os
a = sys.argv[1:]
specs = [x for x in a if ':' in x and not x.startswith('--')]
mn = float(a[a.index('--min') + 1]) if '--min' in a else 0.15
which = a[a.index('--key') + 1] if '--key' in a else 'self'
def norm(s):
    s = re.sub(r'\(.*', '', s) if not s.startswith('<') else s       # C++ args
    s = re.sub(r'^\S+ ', '', s) if re.match(r'^\S+ \S', s) and not s.startswith('<') else s   # return type
    s = s.replace('<rusolve_ekk::', '<')
    m = re.match(r'<([^<>]*)>::(\w+)', s)                        # Rust <path::Class>::method
    if m: cls, meth = m.group(1).split('::')[-1], m.group(2)
    else:
        parts = [p for p in re.split(r'::', s) if p]
        meth = parts[-1] if parts else s
        cls = parts[-2] if len(parts) > 1 else ''
        meth = re.sub(r'\W.*', '', meth)
    return (cls + '.' + meth).lower().replace('_', '')
agg = {'port': collections.defaultdict(float), 'cpp': collections.defaultdict(float)}
for arm in agg:
    for sp in specs:
        i, n = sp.split(':')
        d = json.load(open(f"{os.environ.get('NODEPROF_DIR', 'out')}/{arm}-{i}-{n}.json"))
        tot = d['n']
        rows = d['self'] if which == 'self' else d['incl']
        by = collections.defaultdict(float)
        for name, v in rows: by[norm(name)] += v / tot
        for k, v in by.items(): agg[arm][k] += 100 * v / len(specs)
keys = set(agg['port']) | set(agg['cpp'])
rows = sorted(keys, key=lambda k: -(agg['port'][k] - agg['cpp'][k]))
print(f'mean {which} share over {len(specs)} instances (percent of run); top port-minus-cpp')
print(f'{"function":60s} {"port%":>7s} {"cpp%":>7s} {"diff":>6s}')
for k in rows:
    p, c = agg['port'][k], agg['cpp'][k]
    if abs(p - c) >= mn: print(f'{k[:60]:60s} {p:7.2f} {c:7.2f} {p-c:6.2f}')
