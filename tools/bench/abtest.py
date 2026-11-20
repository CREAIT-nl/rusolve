"""Parallel MIP A/B harness: does an option or a code change make the solver faster?

Why this exists: wall-time A/Bs on a handful of MIPLIB instances cannot see effects
under ~12 % and cost minutes per arm. This harness makes an A/B cheap enough to
run per decision:
  * the metric is **instructions retired** (macOS `/usr/bin/time -l`), which
    barely moves with machine load, so many runs go in parallel on one machine
    (`--metric time` reports solver time instead);
  * instances come from a calibrated **pool** (HiGHS solves each in about
    0.5-30 s), so no CPU is burnt on censored time-limit runs;
  * `probe` finds the instances a flag actually changes (node-limited runs
    compared on status/objective/bound/nodes/LP iterations), so a flag's A/B
    runs only on its **affected set**; unaffected instances are exactly 1.00x;
  * `run --race` stops measuring an arm once its bootstrap CI has decided.

Arms are `name=spec`; spec is `kind[+key=value...]`:
  port     the rusolve-ekk HiGHS port, CLI default profile (HighsOptions::recommended)
  port115  the port with RUSOLVE_PROFILE=highs-1.15.1 (every rusolve-specific option off)
  highspy  upstream HiGHS through highspy (tools/bench/abtest_highs.py), run by the
           `python` that runs this script, or by $ABTEST_HIGHSPY_PYTHON (e.g. a venv
           with another highspy version)
  hbin:<path>  a HiGHS `highs` executable built from source
  e.g. `conf=port+mip_conflict_resolution=true`.

Instances: pool files list paths under ~/.cache/rusolve-bench/miplib-benchmark/
(MIPLIB 2017 `.mps.gz` files, see tools/bench/identity.sh); set
$RUSOLVE_MIPLIB_DIR to read them from another directory. `report` and `pool`
check claimed optima against the MIPLIB `miplib2017.solu` file in that directory.

Subcommands:
  calibrate  run one arm over instance files, write per-instance times
  pool       select pool instances from a calibrate CSV
  run        full solves of arms x pool x seeds  (resumable: rows already in --out are skipped)
  probe      node-limited solves; then `affected` lists, per arm, the instances it changes
  report     paired comparison of every arm against --base

Typical:
  python tools/bench/abtest.py run --pool tools/bench/pools/t2.txt --seeds 0,1 \
      --arm base=port --arm conf=port+mip_conflict_resolution=true --out ab.csv
  python tools/bench/abtest.py report ab.csv --base base
"""

import argparse
import csv
import math
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
HIGHSPY_PYTHON = os.environ.get("ABTEST_HIGHSPY_PYTHON", sys.executable)
HIGHS_RUNNER = os.path.join(HERE, "abtest_highs.py")
POOL_DIR_DEFAULT = "~/.cache/rusolve-bench/miplib-benchmark"
MIPLIB_DIR = os.path.expanduser(os.environ.get("RUSOLVE_MIPLIB_DIR", POOL_DIR_DEFAULT))
SOLU_DEFAULT = os.path.join(MIPLIB_DIR, "miplib2017.solu")
FIELDS = ["arm", "inst", "seed", "status", "obj", "bound", "nodes", "lp_iters", "time_s",
          "wall", "instr", "rss_mb", "rc"]
# Instructions per second of the port on a recent Apple-silicon core: ~21e9/s.
# The shift of the shifted geomean is about 1 s worth of work.
INSTR_SHIFT = 2e10
OPT_DIR = tempfile.mkdtemp(prefix="abtest-opts-")
# `highs` app MIP solving report (HighsMipSolver.cpp, v1.15/1.16)
REPORT_RE = {k: re.compile(r"^\s+" + v + r"\s+(\S+)", re.M) for k, v in
             dict(status="Status", obj="Primal bound", bound="Dual bound", nodes="Nodes",
                  lp_iters="LP iterations", time_s="Timing").items()}
REPORT_STATUS = {"Optimal": "Optimal", "Time": "TimeLimit", "Node": "NodeLimit",
                 "Infeasible": "Infeasible"}


def inst_name(path):
    b = os.path.basename(path)
    for ext in (".gz", ".mps", ".lp"):
        if b.endswith(ext):
            b = b[: -len(ext)]
    return b


def parse_arm(text):
    name, _, spec = text.partition("=")
    parts = spec.split("+")
    kind, opts = parts[0], parts[1:]
    if kind not in ("port", "port115", "highspy") and not kind.startswith("hbin:"):
        raise SystemExit(f"arm {text}: unknown kind {kind}")
    return dict(name=name, kind=kind, opts=opts)


def command(arm, binary, path, tl, seed, node_limit):
    env = dict(os.environ)
    if arm["kind"] in ("port", "port115"):
        env["RUSOLVE_MIP_HIGHS"] = "1"
        if arm["kind"] == "port115":
            env["RUSOLVE_PROFILE"] = "highs-1.15.1"
        else:
            env.pop("RUSOLVE_PROFILE", None)
        cmd = [binary, "mip", path, "--csv", "--time-limit", str(tl), "--seed", str(seed)]
        if node_limit is not None:
            cmd += ["--node-limit", str(node_limit)]
        for o in arm["opts"]:
            cmd += ["--opt", o]
    elif arm["kind"].startswith("hbin:"):
        # a HiGHS `highs` executable built from source (e.g. a commit between two tags)
        opts = ["threads=1", "parallel=off", f"time_limit={tl}", f"random_seed={seed}"]
        if node_limit is not None:
            opts.append(f"mip_max_nodes={node_limit}")
        opts += arm["opts"]
        text = "".join(o.replace("=", " = ", 1) + "\n" for o in opts)
        optfile = os.path.join(OPT_DIR, f"{abs(hash(text)):x}.opt")
        if not os.path.exists(optfile):
            with open(optfile, "w") as fh:
                fh.write(text)
        cmd = [arm["kind"][5:], "--model_file", path, "--options_file", optfile]
    else:
        cmd = [HIGHSPY_PYTHON, "-I", HIGHS_RUNNER, path, str(tl), str(seed)]
        if node_limit is not None:
            cmd.append(str(node_limit))
        cmd += arm["opts"]
    return ["/usr/bin/time", "-l"] + cmd, env


TIME_RE = {"instr": re.compile(r"(\d+)\s+instructions retired"),
           "rss": re.compile(r"(\d+)\s+maximum resident set size"),
           "wall": re.compile(r"([\d.]+) real")}


def solve(arm, binary, path, tl, seed, node_limit):
    cmd, env = command(arm, binary, path, tl, seed, node_limit)
    t0 = time.perf_counter()
    # own process group, so a timeout kills the solver under /usr/bin/time too
    # cwd: the HiGHS app writes Highs.log into its working directory
    p = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                         start_new_session=True, cwd=OPT_DIR)
    try:
        out, err = p.communicate(timeout=tl * 2 + 60)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        out, err = p.communicate()
        rc = -9
    wall = time.perf_counter() - t0
    row = dict(arm=arm["name"], inst=inst_name(path), seed=seed, status="Error", obj="nan",
               bound="nan", nodes="", lp_iters="", time_s="", wall=f"{wall:.3f}", instr="",
               rss_mb="", rc=rc)
    lines = out.strip().splitlines()
    hdr = next((k for k, ln in enumerate(lines) if ln.startswith("status,")), None)
    if hdr is not None and hdr + 1 < len(lines):
        f = lines[hdr + 1].split(",", 6)
        row.update(status=f[0], obj=f[1], bound=f[2], nodes=f[3], lp_iters=f[4], time_s=f[5])
    elif "Solving report" in out:
        rep = out[out.index("Solving report"):]
        for k, rx in REPORT_RE.items():
            m = rx.search(rep)
            if m:
                row[k] = m.group(1)
        row["status"] = REPORT_STATUS.get(row["status"], row["status"])
    elif rc == -9:
        row["status"] = "Killed"
    m = TIME_RE["instr"].search(err)
    if m:
        row["instr"] = m.group(1)
    m = TIME_RE["rss"].search(err)
    if m:
        row["rss_mb"] = f"{int(m.group(1)) / 2**20:.0f}"
    return row


def read_rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh))


def read_pool(path):
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                if line.startswith(POOL_DIR_DEFAULT):
                    line = MIPLIB_DIR + line[len(POOL_DIR_DEFAULT):]
                out.append(os.path.expanduser(line))
    return out


def execute(arms, paths, seeds, tl, node_limit, jobs, out, binary):
    done = {(r["arm"], r["inst"], str(r["seed"])) for r in read_rows(out)}
    work = []
    for s in seeds:
        block = []
        for p in paths:
            order = arms[:]
            random.Random(f"{inst_name(p)}:{s}").shuffle(order)  # interleave arms per pair
            block += [(a, p, s) for a in order if (a["name"], inst_name(p), str(s)) not in done
                      and (a.get("insts") is None or inst_name(p) in a["insts"])]
        work += block
    # Longest instances first would need a calibrate CSV; seed order keeps rounds intact.
    new = not os.path.exists(out)
    lock = threading.Lock()
    fh = open(out, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
    n, t0 = 0, time.time()
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = [ex.submit(solve, a, binary, p, tl, s, node_limit) for a, p, s in work]
        for fut in as_completed(futs):
            r = fut.result()
            with lock:
                w.writerow(r)
                fh.flush()
                n += 1
                if n % 25 == 0 or n == len(work):
                    print(f"  {n}/{len(work)} runs, {time.time() - t0:.0f}s", file=sys.stderr)
    fh.close()


def staged_binary(binary, tmp):
    """Copy the binary so a concurrent cargo build cannot swap it mid-run."""
    dst = os.path.join(tmp, "rusolve-cli")
    shutil.copy2(binary, dst)
    return dst


# ---------------------------------------------------------------- analysis

def read_solu(path):
    ref = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 3 and f[0] in ("=opt=", "=best="):
                    ref[f[1]] = (f[0], float(f[2]))
                elif len(f) >= 2 and f[0] == "=inf=":
                    ref[f[1]] = ("=inf=", None)
    return ref


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def wrong(row, ref, gap=1e-4):
    """Optimal claimed but objective inconsistent with the published value."""
    r = ref.get(row["inst"])
    if row["status"] == "Infeasible":
        return r is not None and r[0] != "=inf="
    if row["status"] != "Optimal" or r is None or r[1] is None:
        return False
    obj, pub = fnum(row["obj"]), r[1]
    tol = max(gap * max(abs(obj), abs(pub)), 1e-6) * 1.01 + 1e-6
    if r[0] == "=opt=":
        return abs(obj - pub) > tol
    return False  # =best= only bounds one side; a better value is not wrong


def solved(row):
    return row["status"] in ("Optimal", "Infeasible")


def metric(row, kind):
    if kind == "instr":
        v = fnum(row["instr"])
        return v, INSTR_SHIFT
    v = fnum(row["time_s"]) if row["time_s"] else fnum(row["wall"])
    return v, 1.0


def bootstrap_ci(pairs, B=2000, seed=0):
    """pairs: {inst: [log ratios]}; resample instances (seeds of one instance are correlated)."""
    insts = list(pairs)
    if not insts:
        return float("nan"), float("nan"), float("nan")

    def stat(keys):
        vals = [v for k in keys for v in pairs[k]]
        return sum(vals) / len(vals)

    point = stat(insts)
    rng = random.Random(seed)
    boots = sorted(stat([rng.choice(insts) for _ in insts]) for _ in range(B))
    lo, hi = boots[int(0.025 * B)], boots[int(0.975 * B) - 1]
    return math.exp(point), math.exp(lo), math.exp(hi)


def compare(rows, base, arm, kind, insts=None):
    by = {}
    for r in rows:
        by[(r["arm"], r["inst"], r["seed"])] = r
    pairs, solved_b, solved_a, n = {}, 0, 0, 0
    for (a, i, s), ra in by.items():
        if a != arm or (insts is not None and i not in insts):
            continue
        rb = by.get((base, i, s))
        if rb is None:
            continue
        n += 1
        solved_a += solved(ra)
        solved_b += solved(rb)
        va, sh = metric(ra, kind)
        vb, _ = metric(rb, kind)
        if math.isnan(va) or math.isnan(vb):
            continue
        # unsolved runs keep their (limit-truncated) work: a lower bound on the true cost
        pairs.setdefault(i, []).append(math.log((va + sh) / (vb + sh)))
    g, lo, hi = bootstrap_ci(pairs)
    return dict(n=n, insts=len(pairs), ratio=g, lo=lo, hi=hi, solved_a=solved_a, solved_b=solved_b)


def report(args):
    rows = read_rows(args.csv)
    ref = read_solu(args.solu)
    arms = sorted({r["arm"] for r in rows})
    bad = [r for r in rows if wrong(r, ref)]
    errs = [r for r in rows if r["status"] in ("Error", "Killed")]
    insts = None
    if args.insts:
        insts = set(read_pool(args.insts)) if os.path.exists(args.insts) else set(args.insts.split(","))
        insts = {inst_name(i) for i in insts}
    print(f"{args.csv}: {len(rows)} rows, metric={args.metric}, base={args.base}")
    print(f"{'arm':28s} {'pairs':>5s} {'inst':>4s} {'ratio':>6s} {'95% CI':>15s} {'solved':>11s}  verdict")
    for a in arms:
        if a == args.base:
            continue
        c = compare(rows, args.base, a, args.metric, insts)
        if c["n"] == 0:
            continue
        v = decide(c, args.band)
        print(f"{a:28s} {c['n']:5d} {c['insts']:4d} {c['ratio']:6.3f} "
              f"[{c['lo']:.3f},{c['hi']:.3f}] {c['solved_a']:4d} vs {c['solved_b']:<4d}  {v}")
    for r in bad:
        print(f"WRONG {r['arm']} {r['inst']} s{r['seed']} {r['status']} obj={r['obj']} "
              f"ref={ref.get(r['inst'])}")
    for r in errs:
        print(f"ERR   {r['arm']} {r['inst']} s{r['seed']} {r['status']} rc={r['rc']}")
    return 1 if bad else 0


def decide(c, band):
    if math.isnan(c["ratio"]):
        return "n/a"
    if c["insts"] < 3:  # the bootstrap resamples instances: <3 gives a degenerate CI
        return "too-few"
    if c["hi"] < 1.0:
        return "WIN"
    if c["lo"] > 1.0:
        return "LOSS"
    if c["lo"] >= 1 - band and c["hi"] <= 1 + band:
        return "neutral"
    return "undecided"


def affected(args):
    rows = read_rows(args.csv)
    by = {(r["arm"], r["inst"], r["seed"]): r for r in rows}
    key = ("status", "obj", "bound", "nodes", "lp_iters")
    arms = sorted({r["arm"] for r in rows} - {args.base})
    insts = sorted({r["inst"] for r in rows})
    out = {}
    for a in arms:
        hit = []
        for i in insts:
            for (aa, ii, s), ra in by.items():
                if aa != a or ii != i:
                    continue
                rb = by.get((args.base, i, s))
                if rb and any(ra[k] != rb[k] for k in key):
                    hit.append(i)
                    break
        out[a] = hit
        print(f"{a:34s} {len(hit):3d}/{len(insts)}  {' '.join(hit)}")
    if args.write_dir:
        os.makedirs(args.write_dir, exist_ok=True)
        for a, hit in out.items():
            with open(os.path.join(args.write_dir, f"{a}.txt"), "w") as fh:
                fh.write("\n".join(hit) + ("\n" if hit else ""))


def pool(args):
    rows = [r for r in read_rows(args.csv) if r["arm"] == args.arm]
    ref = read_solu(args.solu)
    picked = []
    for r in sorted(rows, key=lambda r: fnum(r["time_s"])):
        t = fnum(r["time_s"])
        if solved(r) and not wrong(r, ref) and args.tmin <= t <= args.tmax:
            picked.append((r["inst"], t, r["nodes"]))
    paths = {inst_name(p): p for d in args.dirs for p in instance_files(d)}
    with open(args.write, "w") as fh:
        fh.write(f"# {len(picked)} instances: {args.arm} solve time in [{args.tmin},{args.tmax}] s, "
                 f"threads=1, seed 0 (from {os.path.basename(args.csv)})\n")
        for i, t, nodes in picked:
            fh.write(f"{paths[i]}  # {t:.2f}s {nodes} nodes\n")
    print(f"{len(picked)} instances -> {args.write}")


def instance_files(d):
    d = os.path.expanduser(d)
    seen = {}
    for f in sorted(os.listdir(d)):
        if f.endswith((".mps", ".mps.gz")):
            n = inst_name(f)
            # prefer plain .mps over .mps.gz for the same instance
            if n not in seen or f.endswith(".mps"):
                seen[n] = os.path.join(d, f)
    return list(seen.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def run_args(p, default_tl):
        p.add_argument("--arm", action="append", required=True)
        p.add_argument("--pool", help="pool file (one path per line)")
        p.add_argument("--dir", action="append", default=[], help="instance directory (calibrate)")
        p.add_argument("--insts", help="restrict to these instance names (file or comma list)")
        p.add_argument("--seeds", default="0")
        p.add_argument("--time-limit", type=float, default=default_tl)
        p.add_argument("--jobs", type=int, default=12)
        p.add_argument("--binary", default=os.path.join(REPO, "target", "release", "rusolve-cli"))
        p.add_argument("--out", required=True)

    for name, tl in (("calibrate", 60.0), ("run", 60.0), ("probe", 60.0)):
        p = sub.add_parser(name)
        run_args(p, tl)
        if name == "probe":
            p.add_argument("--node-limit", type=int, default=50)
        if name == "run":
            p.add_argument("--affected-dir",
                           help="run each non-base arm only on the instances in <dir>/<arm>.txt")
            p.add_argument("--race", action="store_true",
                           help="after each seed round, drop arms whose CI has decided")
            p.add_argument("--base", default="base")
            p.add_argument("--band", type=float, default=0.03)
            p.add_argument("--metric", default="instr", choices=["instr", "time"])
    p = sub.add_parser("report")
    p.add_argument("csv")
    p.add_argument("--base", default="base")
    p.add_argument("--metric", default="instr", choices=["instr", "time"])
    p.add_argument("--band", type=float, default=0.03, help="CI inside 1±band counts as neutral")
    p.add_argument("--insts")
    p.add_argument("--solu", default=SOLU_DEFAULT)
    p = sub.add_parser("affected")
    p.add_argument("csv")
    p.add_argument("--base", default="base")
    p.add_argument("--write-dir")
    p = sub.add_parser("pool")
    p.add_argument("csv")
    p.add_argument("--arm", default="highspy")
    p.add_argument("--tmin", type=float, default=0.5)
    p.add_argument("--tmax", type=float, default=30.0)
    p.add_argument("--dir", dest="dirs", action="append", required=True)
    p.add_argument("--solu", default=SOLU_DEFAULT)
    p.add_argument("--write", required=True)
    args = ap.parse_args()

    if args.cmd == "report":
        sys.exit(report(args))
    if args.cmd == "affected":
        return affected(args)
    if args.cmd == "pool":
        return pool(args)

    arms = [parse_arm(a) for a in args.arm]
    if args.pool:
        paths = read_pool(args.pool)
    else:
        uniq = {}
        for d in args.dir:  # earlier directories win for a duplicated instance name
            for p in instance_files(d):
                uniq.setdefault(inst_name(p), p)
        paths = list(uniq.values())
    if args.insts:
        keep = set(read_pool(args.insts)) if os.path.exists(args.insts) else set(args.insts.split(","))
        keep = {inst_name(k) for k in keep}
        paths = [p for p in paths if inst_name(p) in keep]
    if getattr(args, "affected_dir", None):
        union = set()
        for a in arms:
            if a["name"] == args.base:
                continue
            a["insts"] = set(read_pool(os.path.join(args.affected_dir, f"{a['name']}.txt")))
            union |= a["insts"]
        arms = [a for a in arms if a["name"] == args.base or a["insts"]]
        paths = [p for p in paths if inst_name(p) in union]
    seeds = [int(s) for s in args.seeds.split(",")]
    node_limit = getattr(args, "node_limit", None)
    with tempfile.TemporaryDirectory(prefix="abtest-") as tmp:
        binary = staged_binary(args.binary, tmp) if any(a["kind"].startswith("port") for a in arms) else None
        if args.cmd == "run" and args.race:
            live = arms[:]
            for k, s in enumerate(seeds):
                execute(live, paths, [s], args.time_limit, None, args.jobs, args.out, binary)
                rows = read_rows(args.out)
                keep = []
                for a in live:
                    if a["name"] == args.base:
                        keep.append(a)
                        continue
                    c = compare(rows, args.base, a["name"], args.metric)
                    v = decide(c, args.band)
                    print(f"  seed {s}: {a['name']} {c['ratio']:.3f} [{c['lo']:.3f},{c['hi']:.3f}] {v}",
                          file=sys.stderr)
                    # one seed round never decides: a single seed is one tree draw per instance
                    if k == 0 or v not in ("WIN", "LOSS", "neutral"):
                        keep.append(a)
                live = keep
                if len(live) <= 1:
                    break
        else:
            execute(arms, paths, seeds, args.time_limit, node_limit, args.jobs, args.out, binary)


if __name__ == "__main__":
    main()
