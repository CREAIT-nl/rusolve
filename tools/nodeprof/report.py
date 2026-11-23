#!/usr/bin/env python3
"""Resolve a nodeprof dump to per-function and per-subsystem cycle shares.
usage: report.py prof.bin exe [--rules rules.txt] [--top N] [--json out.json]
Symbols come from `nm -C -n` on exe (PIE load base from the dumped maps). Subsystem of a
sample = the rule matched by its innermost matching frame (file order breaks ties within a frame)."""
import bisect, collections, json, re, struct, subprocess, sys

def load(path):
    raw = open(path, 'rb').read()
    magic, depth, n, lost = struct.unpack_from('<4I', raw, 0)
    assert magic == 0x4e50524f
    samples = struct.unpack_from('<%dQ' % (depth * n), raw, 16)
    maps = raw[16 + 8 * depth * n:].decode(errors='replace').splitlines()
    return depth, n, lost, [samples[i * depth:(i + 1) * depth] for i in range(n)], maps

def symtab(exe):
    out = subprocess.run(['nm', '-C', '-n', '--defined-only', exe], capture_output=True, text=True).stdout
    addrs, names = [], []
    for l in out.splitlines():
        p = l.split(' ', 2)
        if len(p) == 3 and p[1] in 'tTwW':
            addrs.append(int(p[0], 16)); names.append(p[2])
    return addrs, names

def main():
    a = sys.argv[1:]
    prof, exe = a[0], a[1]
    rules_f = a[a.index('--rules') + 1] if '--rules' in a else None
    top = int(a[a.index('--top') + 1]) if '--top' in a else 25
    jout = a[a.index('--json') + 1] if '--json' in a else None
    only = a[a.index('--only') + 1] if '--only' in a else None   # restrict self/incl tables to one subsystem (prefix)
    depth, n, lost, samples, maps = load(prof)
    # modules: (lo, hi, offset, path) from the dumped maps; the main exe is resolved with
    # nm on `exe`, other mapped files (libc ...) with nm -D, else labelled [lib].
    mods, bases = [], {}
    for l in maps:
        p = l.split()
        if len(p) >= 6 and p[5].startswith('/'):
            lo, hi = (int(x, 16) for x in p[0].split('-'))
            bases[p[5]] = min(bases.get(p[5], 1 << 62), lo - int(p[2], 16))
            if 'x' in p[1]: mods.append((lo, hi, int(p[2], 16), p[5]))
    exename = exe.split('/')[-1]
    tabs = {}
    def tab(path):
        if path not in tabs:
            if path.endswith(exename): tabs[path] = symtab(exe)
            else:
                out = subprocess.run(['nm', '-D', '-n', '--defined-only', path], capture_output=True, text=True).stdout
                ad, nm = [], []
                for l in out.splitlines():
                    q = l.split(' ', 2)
                    if len(q) == 3 and q[1] in 'tTwWiI': ad.append(int(q[0], 16)); nm.append(q[2])
                tabs[path] = (ad, nm)
        return tabs[path]
    pie = subprocess.run(['file', '-L', exe], capture_output=True, text=True).stdout
    exe_pie = 'pie' in pie or 'shared object' in pie
    cache = {}
    def sym(ip):
        if ip in cache: return cache[ip]
        s = '?'
        for lo, hi, off, path in mods:
            if lo <= ip < hi:
                is_exe = path.endswith(exename)
                r = ip - bases[path]
                if is_exe and not exe_pie: r = ip
                addrs, names = tab(path)
                i = bisect.bisect_right(addrs, r) - 1
                if i >= 0 and r - addrs[i] < (0x400 if not is_exe else (1 << 22)): s = names[i]
                else: s = '[' + path.split('/')[-1] + ']'
                if not is_exe and s.startswith('['): pass
                break
        s = re.sub(r'\s*\[clone .*?\]', '', s)
        s = re.sub(r'<.*', '', s) if s.count('<') > 3 else s
        cache[ip] = s
        return s
    rules = []
    if rules_f:
        for l in open(rules_f):
            l = l.strip()
            if l and not l.startswith('#'):
                name, pat = l.split(None, 1); rules.append((name, re.compile(pat.lower().replace('_', ''))))
    selfc, incl, sub = collections.Counter(), collections.Counter(), collections.Counter()
    libc_by = collections.Counter()
    exe_rng = [(lo, hi) for lo, hi, off, path in mods if path.endswith(exename)]
    in_exe = lambda ip: any(lo <= ip < hi for lo, hi in exe_rng)
    kept = 0
    for s in samples:
        s = list(s)
        if in_exe(s[0]) or not in_exe(s[1]): del s[1]
        fr = [sym(ip) for ip in s if ip]
        tag = 'other'
        norm = [f.lower().replace('_', '') for f in fr]
        done = False
        for f in norm:
            for name, pat in rules:
                if pat.search(f): tag = name; done = True; break
            if done: break
        sub[tag] += 1
        if '.' in tag: sub[tag.split('.')[0] + '.*'] += 1
        if only and not tag.startswith(only): continue
        kept += 1
        selfc[fr[0]] += 1
        if fr[0].startswith('['):
            c = next((f for f in fr if not f.startswith('[')), '?'); libc_by[c] += 1
        for f in set(fr): incl[f] += 1
    if only: n_all = n; n = kept; print(f'(restricted to {only}: {kept} of {n_all} samples; shares below are of the restricted set)')
    print(f'samples {n} lost {lost}')
    print('-- subsystem (innermost matching frame)')
    for k, v in sorted(sub.items()):
        if v: print(f'{100*v/n:6.2f}%  {k}')
    print('-- self'); 
    for k, v in selfc.most_common(top): print(f'{100*v/n:6.2f}%  {k[:150]}')
    print('-- self in unresolved library code (libc memcpy/malloc/...), by nearest resolved caller')
    for k, v in libc_by.most_common(top): print(f'{100*v/n:6.2f}%  {k[:150]}')
    print('-- inclusive')
    for k, v in incl.most_common(top): print(f'{100*v/n:6.2f}%  {k[:150]}')
    if jout:
        json.dump({'n': n, 'sub': sub, 'self': selfc.most_common(200), 'incl': incl.most_common(200)}, open(jout, 'w'))
main()
