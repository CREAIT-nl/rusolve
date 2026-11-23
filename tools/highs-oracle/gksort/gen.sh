#!/bin/sh
# Regenerate the GK_MKQSORT reference outputs the independent Rust sort is tested against.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../../.." && pwd)
highs=${HIGHS_SRC:-$root/reference/HiGHS}
[ -d "$highs/extern/metis/GKlib" ] || { echo "gen.sh: set HIGHS_SRC to a HiGHS checkout (extern/metis/GKlib not found under $highs)" >&2; exit 1; }
out=$root/crates/rusolve-ekk/tests/fixtures/gksort
mkdir -p "$out"
tmp=$(mktemp -d)
cc -O1 -std=c99 -I"$highs/extern/metis/GKlib" "$here/gksort_ref.c" -o "$tmp/gksort_ref"
"$tmp/gksort_ref" > "$out/gksort.txt"
rm -rf "$tmp"
wc -lc "$out/gksort.txt"
