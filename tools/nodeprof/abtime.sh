#!/bin/bash
# usage: tools/nodeprof/abtime.sh <instance> <nodes> <reps> <binA> <binB> ; interleaved A,B,A,B user-time
# runs of two rusolve-cli binaries (or "cpp" for a HiGHS build, $HIGHS_BIN, default
# target/highs16/bin/highs). Prints user seconds per run and the ratio. Run from the checkout root
# (Linux, GNU time); reads $NODEPROF_DIR/<instance>.mps and $NODEPROF_DIR/faithful.opt (default out/).
D=${NODEPROF_DIR:-out}
HIGHS_BIN=${HIGHS_BIN:-target/highs16/bin/highs}
inst=$1 nodes=$2 reps=$3 A=$4 B=$5
run() { # bin -> prints "user nodes lpit obj"
  if [ "$1" = cpp ]; then
    printf "threads = 1\nrandom_seed = 0\nmip_max_nodes = $nodes\n" > $D/ab-$inst.opt
    /usr/bin/time -f "user=%U" -o $D/ab.t $HIGHS_BIN $D/$inst.mps --options_file $D/ab-$inst.opt > $D/ab.o 2>&1
    echo "$(grep user= $D/ab.t | sed "s/user=//") $(awk '/^  Nodes/{n=$2} /LP iterations/{i=$3} END{print n, i}' $D/ab.o)"
  else
    RUSOLVE_MIP_HIGHS=1 /usr/bin/time -f "user=%U" -o $D/ab.t $1 mip $D/$inst.mps --csv --options-file $D/faithful.opt --node-limit $nodes > $D/ab.o 2>/dev/null
    echo "$(sed 's/user=//' $D/ab.t) $(tail -1 $D/ab.o | cut -d, -f1-5 | awk -F, '{print $4, $5, $2}')"
  fi
}
ulimit -v 6000000
for r in $(seq $reps); do
  echo "A $(run $A)"; echo "B $(run $B)"
done | awk -v inst=$inst '{ if ($1=="A") {a[++na]=$2; ra=$3" "$4" "$5} else {b[++nb]=$2; rb=$3" "$4" "$5} }
 END { sa=0; sb=0; for(i=1;i<=na;i++){sa+=a[i]; sb+=b[i]; printf "%s rep%d A %.2f B %.2f ratio B/A %.3f\n", inst,i,a[i],b[i],b[i]/a[i]}
  printf "%s MEAN A %.2f B %.2f B/A %.3f   A:[%s] B:[%s]\n", inst, sa/na, sb/nb, sb/sa, ra, rb }'
