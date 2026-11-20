#!/bin/bash
# usage: tools/nodeprof/prof.sh <arm: port|cpp|portplain|cppplain> <instance> <node_limit> [time_limit]
# Runs one solve under the nodeprof sampler (port|cpp) or plain with user-time (portplain|cppplain),
# to compare where the Rust port and the HiGHS C++ build spend their time (Linux, GNU time).
# Run from the checkout root. Setup:
#   cc -O2 -shared -fPIC -o tools/nodeprof/sampler.so tools/nodeprof/sampler.c
#   $NODEPROF_DIR/<instance>.mps      the model (default directory: out/)
#   $NODEPROF_DIR/faithful.opt        options file passed to the port (may be empty)
#   PORT_BIN   rusolve-cli built with frame pointers (default target-fp/release/rusolve-cli:
#              RUSTFLAGS="-C force-frame-pointers=yes" cargo build --release -p rusolve-cli --target-dir target-fp)
#   HIGHS_BIN  a HiGHS `highs` executable built with -fno-omit-frame-pointer (default target/highs16/bin/highs)
# Output: $NODEPROF_DIR/<arm>-<inst>-<nodes>.{prof,res,err}.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
arm=$1 inst=$2 nodes=$3 tl=${4:-600}
D=${NODEPROF_DIR:-out}; tag=$arm-$inst-$nodes
PORT_BIN=${PORT_BIN:-target-fp/release/rusolve-cli}
HIGHS_BIN=${HIGHS_BIN:-target/highs16/bin/highs}
ulimit -v 6000000
printf "threads = 1\nrandom_seed = 0\nmip_max_nodes = $nodes\ntime_limit = $tl\n" > $D/$tag.opt
case $arm in
 port)
  RUSOLVE_MIP_HIGHS=1 /usr/bin/time -f "user=%U sys=%S maxrss=%M" env NODEPROF_OUT=$D/$tag.prof LD_PRELOAD=$HERE/sampler.so $PORT_BIN mip $D/$inst.mps --csv \
   --options-file $D/faithful.opt --node-limit $nodes --time-limit $tl > $D/$tag.res 2> $D/$tag.err ;;
 portplain)
  RUSOLVE_MIP_HIGHS=1 /usr/bin/time -f "user=%U sys=%S maxrss=%M" ${BIN:-target/release/rusolve-cli} mip $D/$inst.mps --csv \
   --options-file $D/faithful.opt --node-limit $nodes --time-limit $tl > $D/$tag.res 2> $D/$tag.err ;;
 cppplain)
  /usr/bin/time -f "user=%U sys=%S maxrss=%M" $HIGHS_BIN $D/$inst.mps \
   --options_file $D/port-$inst-$nodes.opt > $D/$tag.res 2> $D/$tag.err ;;
 cpp)
  /usr/bin/time -f "user=%U sys=%S maxrss=%M" env NODEPROF_OUT=$D/$tag.prof LD_PRELOAD=$HERE/sampler.so $HIGHS_BIN $D/$inst.mps \
   --options_file $D/$tag.opt > $D/$tag.res 2> $D/$tag.err ;;
esac
tail -1 $D/$tag.err
