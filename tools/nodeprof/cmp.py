#!/usr/bin/env python3
"""Side-by-side port vs C++ subsystem seconds for one instance.
usage: tools/nodeprof/cmp.py <inst> <nodes>   (run from the checkout root, after prof.sh port and prof.sh cpp)
Seconds = sample share x user CPU time of the run (from /usr/bin/time in the .err file).
Reads $NODEPROF_DIR (default out/); binaries from $PORT_BIN and $HIGHS_BIN as in prof.sh."""
import json, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.environ.get('NODEPROF_DIR', 'out')
inst, nodes = sys.argv[1:3]
def run(arm, exe):
    tag = f'{D}/{arm}-{inst}-{nodes}'
    subprocess.run([sys.executable, os.path.join(HERE, 'report.py'), tag + '.prof', exe, '--rules', os.path.join(HERE, 'rules.txt'), '--top', '0', '--json', tag + '.json'],
                   capture_output=True)
    d = json.load(open(tag + '.json'))
    err = open(tag + '.err').read()
    user = float(re.search(r'user=([\d.]+)', err).group(1))
    res = open(tag + '.res').read()
    if arm == 'port':
        f = res.strip().splitlines()[-1].split(',')
        nd, it = int(f[3]), int(f[4])
    else:
        nd = int(re.search(r'^\s*Nodes\s+(\d+)', res, re.M).group(1))
        it = int(re.search(r'LP iterations\s+(\d+)', res).group(1))
    n = d['n']
    sec = {k: v / n * user for k, v in d['sub'].items()}
    return user, nd, it, sec
pu, pn, pi, ps = run('port', os.environ.get('PORT_BIN', 'target-fp/release/rusolve-cli'))
cu, cn, ci, cs = run('cpp', os.environ.get('HIGHS_BIN', 'target/highs16/bin/highs'))
print(f'{inst} nodes {nodes}: port user {pu:.1f}s nodes {pn} lpit {pi} | cpp user {cu:.1f}s nodes {cn} lpit {ci}')
print(f'{"subsystem":18s} {"port s":>8s} {"cpp s":>8s} {"ratio":>6s} {"port/it ms":>10s} {"cpp/it ms":>10s}')
for k in sorted(set(ps) | set(cs)):
    a, b = ps.get(k, 0), cs.get(k, 0)
    if a + b < 0.05: continue
    print(f'{k:18s} {a:8.2f} {b:8.2f} {a/b if b else float("inf"):6.2f} {1000*a/pi:10.3f} {1000*b/ci:10.3f}')
print(f'{"TOTAL":18s} {pu:8.2f} {cu:8.2f} {pu/cu:6.2f} {1000*pu/pi:10.3f} {1000*cu/ci:10.3f}')
