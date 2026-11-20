#!/bin/zsh
# Parity check for contributors: does my change alter the default search path?
#
# Solves a pool of MIPLIB instances with two rusolve-cli builds (typically one built from main and
# one from your branch) and compares the results. A refactor, a performance change, debug-only code
# or a new default-off option must leave every instance identical.
#
# usage: tools/bench/identity.sh <cli_A> <cli_B> [pool] [--opt k=v ...]
#   pool defaults to tools/bench/pools/t2fast.txt (20 MIPLIB 2017 instances that solve fast).
#   Extra args go to both binaries (e.g. --opt mip_heuristic_run_local_mip=true).
#   IDENTITY_B_OPTS="--opt k=v ..." goes to <cli_B> only: switch a default-on option off on the
#   new binary so it can still be compared with an older reference that lacks the option.
# Compares the first five CSV fields (status, primal, dual, nodes, lp_iters) of
#   RUSOLVE_MIP_HIGHS=1 <cli> mip <inst> --csv --time-limit 120 --seed 0
# runs 4 at a time and prints "IDENTICAL n of m" plus every differing instance.
#
# Instances: <name>.mps.gz files in $RUSOLVE_MIPLIB_DIR (default
# ~/.cache/rusolve-bench/miplib-benchmark), downloaded from MIPLIB 2017, e.g.
#   mkdir -p ~/.cache/rusolve-bench/miplib-benchmark && cd ~/.cache/rusolve-bench/miplib-benchmark
#   for i in $(cat <repo>/tools/bench/pools/t2fast.txt); do
#     curl -fLO https://miplib.zib.de/WebData/instances/$i.mps.gz; done
# A time-limited instance is no witness: compare it under a node limit instead.
set -e
A=$1 B=$2; shift 2
pool=tools/bench/pools/t2fast.txt
if [[ $# -gt 0 && $1 != --* ]]; then pool=$1; shift; fi
extra=("$@")
dir=${RUSOLVE_MIPLIB_DIR:-$HOME/.cache/rusolve-bench/miplib-benchmark}
one() {
  local arm=$1 bin=$2 i=$3
  local -a more=()
  [[ $arm == B ]] && eval "more=(${IDENTITY_B_OPTS:-})"
  echo "$arm,$i,$(RUSOLVE_MIP_HIGHS=1 $bin mip $dir/$i.mps.gz --csv --time-limit 120 --seed 0 "${extra[@]}" "${more[@]}" 2>/dev/null | tail -1 | cut -d, -f1-5)"
}
out=$(mktemp)
for i in $(cat $pool); do
  one A $A $i >> $out & one B $B $i >> $out &
  while [[ $(jobs -r | wc -l) -ge 4 ]]; do sleep 1; done
done
wait
python3 - "$out" <<'PY'
import collections, sys
d = collections.defaultdict(dict)
for l in open(sys.argv[1]):
    arm, i, rest = l.strip().split(',', 2); d[i][arm] = rest
bad = [i for i in d if d[i].get('A') != d[i].get('B')]
print('IDENTICAL', len(d) - len(bad), 'of', len(d))
for i in bad: print('  DIFF', i, 'A=', d[i].get('A'), 'B=', d[i].get('B'))
PY
rm -f $out
