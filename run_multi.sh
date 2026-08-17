#!/bin/sh
# run_multi.sh N MU [OUTDIR]  -- multigraph collision search, one edge bucket
# at a time.  Collisions occur only within a fixed edge count, so bucketing is
# lossless and keeps peak memory to a single bucket instead of the whole range.
#
#   sh run_multi.sh 6 4
#
# Needs nauty's geng and multig on the PATH, and an mg_collisions built from a
# source that supports -m (rebuild if in doubt: make clean && make).
set -e
N=${1:?usage: run_multi.sh N MU [OUTDIR]}
MU=${2:?usage: run_multi.sh N MU [OUTDIR]}
OUT=${3:-res_${N}_${MU}}

# ---- fail early, and clearly, on a stale binary -------------------------
if ! ./mg_collisions --help 2>&1 | grep -q '\-m EDGES'; then
  echo "ERROR: ./mg_collisions does not support -m." >&2
  echo "       Rebuild from the current source:  make clean && make" >&2
  exit 1
fi

mkdir -p "$OUT"
MAXM=$(( N*(N-1)/2*MU ))
if [ -s "$OUT/all.mg" ]; then
  echo "reusing $OUT/all.mg"
else
  echo "generating multigraphs ..."
  geng -q -c "$N" | multig -m"$MU" -T > "$OUT/all.mg"
fi

for m in $(seq $((N-1)) $MAXM); do
  [ -s "$OUT/m$m.txt" ] && continue
  printf "m=%s " "$m"
  ./mg_collisions -q -m "$m" < "$OUT/all.mg" > "$OUT/m$m.txt"
  grep -h "^SUMMARY L " "$OUT/m$m.txt" | sed 's/SUMMARY L //'
done

echo
echo "totals:"
awk '/^SUMMARY/{split($4,a,"=");split($5,b,"=");s[$2]+=b[2];g[$2]+=a[2]}
     END{for(k in s) printf "  %-5s graphs=%d classes=%d\n",k,g[k],s[k]}' "$OUT"/m*.txt
echo
echo "checks that MUST hold in every bucket (any line below is a finding):"
# L and Lhat must induce the same partition -- this is the positive statement
grep -h "L   -classes refine Lhat" "$OUT"/m*.txt | grep -v ": YES" || true
grep -h "Lhat-classes refine L" "$OUT"/m*.txt | grep -v ": YES" || true
grep -h "L-classes    == Lhat" "$OUT"/m*.txt | grep -v ": YES" || true
# the surviving implication of the conjecture
grep -h "L   -classes refine XB" "$OUT"/m*.txt | grep -v ": YES" || true
grep -h "Lhat-classes refine XB" "$OUT"/m*.txt | grep -v ": YES" || true
# Sarmiento: U and XB are equivalent for all graphs, so this is a code check
grep -h "XB  -classes refine U" "$OUT"/m*.txt | grep -v ": YES" || true
grep -h "U   -classes refine XB" "$OUT"/m*.txt | grep -v ": YES" || true
echo "(end of must-hold checks)"

echo
echo "expected to FAIL somewhere -- this is the multigraph counterexample:"
printf "  XB refines L  : %s\n" \
  "$(grep -h 'XB  -classes refine L' "$OUT"/m*.txt | grep -c ': NO') bucket(s) with NO"
printf "  Lhat == XB    : %s\n" \
  "$(grep -h 'Lhat-classes == XB' "$OUT"/m*.txt | grep -c ': NO') bucket(s) with NO"
