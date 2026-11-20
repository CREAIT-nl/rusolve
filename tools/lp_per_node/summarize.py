#!/usr/bin/env python3
"""Summarise stats.phases.lp_iterations.by_source of rusolve mip --format json results.
usage: summarize.py [--depth] result.json...   (markdown tables)"""
import json, sys, re

def groups(by):
    g = {}
    def add(k, n):
        it, so = g.get(k, (0, 0)); g[k] = (it + n['iterations'], so + n['solves'])
    for name, n in by.items():
        m = re.match(r'node_(child|sibling|jump)_', name)
        if m: add('node_' + m.group(1), n); add('node_all', n)
        elif name.startswith('heur_'): add('heur_all', n); add(name, n)
        else: add(name, n)
        add('all', n)
    return g

def load(p):
    d = json.load(open(p)); st = d['stats']; ph = st.get('phases')
    if not ph or 'by_source' not in ph['lp_iterations']: return None
    return d, st, ph, groups(ph['lp_iterations']['by_source'])

def name(p): return re.sub(r'^(.*/)?[a-z0-9]+-(.*)-\d+\.json$', r'\2', p)

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    depth = '--depth' in sys.argv
    print('| instance | nodes | reported total | true sum | root | root sepa | tree sepa | strong br | node child | node sibling | node jump | heuristics | untagged |')
    print('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    rows = []
    for p in args:
        r = load(p)
        if not r: continue
        d, st, ph, g = r; f = lambda k: g.get(k, (0, 0))[0]
        rows.append((p, st, ph, g))
        print(f"| {name(p)} | {st['nodes']} | {ph['lp_iterations']['total']} | {f('all')} | {f('root_lp')} | {f('root_separation')} | {f('tree_separation')} | {f('strong_branching')} | {f('node_child')} | {f('node_sibling')} | {f('node_jump')} | {f('heur_all')} | {f('untagged')} |")
    print()
    print('Per node (iterations / nodes) and per solve:')
    print('| instance | node LPs | node it/node | node it/LP | child it/LP | sibling it/LP | jump it/LP (n) | tree sepa/node | SB/node | heur/node | all/node |')
    print('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for p, st, ph, g in rows:
        n = max(st['nodes'], 1); G = lambda k: g.get(k, (0, 0))
        ni, ns = G('node_all')
        per = lambda k: f"{G(k)[0] / G(k)[1]:.1f}" if G(k)[1] else '-'
        print(f"| {name(p)} | {ns} | {ni / n:.1f} | {ni / max(ns, 1):.1f} | {per('node_child')} | {per('node_sibling')} | {per('node_jump')} ({G('node_jump')[1]}) | {G('tree_separation')[0] / n:.1f} | {G('strong_branching')[0] / n:.1f} | {G('heur_all')[0] / n:.1f} | {G('all')[0] / n:.1f} |")
    if depth:
        print()
        for p, st, ph, g in rows:
            print(f'{name(p)} node re-solves by origin and depth (iterations / solves = it per LP):')
            for o in ('child', 'sibling', 'jump'):
                cells = []
                for b in ('0_5', '6_15', '16_30', '31_plus'):
                    it, so = g.get(f'node_{o}_depth_{b}', (0, 0))
                    cells.append(f'{it / so:.1f} ({so})' if so else '-')
                print(f'  {o:8s} ' + '  '.join(cells))
