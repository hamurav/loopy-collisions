#!/bin/sh
# check_partial.sh OUTDIR  -- run the must-hold checks on whatever buckets exist
OUT=${1:-res_6_4}
echo "buckets finished: $(ls "$OUT"/m*.txt 2>/dev/null | wc -l)"
echo
echo "MUST hold (any line below is a genuine finding):"
for pat in "L   -classes refine Lhat" "Lhat-classes refine L" \
           "L-classes    == Lhat" "L   -classes refine XB" \
           "Lhat-classes refine XB" "XB  -classes refine U" "U   -classes refine XB"; do
  grep -h "$pat" "$OUT"/m*.txt 2>/dev/null | grep -v ": YES"
done
echo "(nothing above = all good)"
echo
echo "running totals so far:"
awk '/^SUMMARY/{split($4,a,"=");split($5,b,"=");split($6,p,"=");
                s[$2]+=b[2];g[$2]+=a[2];q[$2]+=p[2]}
     END{for(k in s) printf "  %-5s graphs=%-9d classes=%-6d pairs=%d\n",k,g[k],s[k],q[k]}' \
    "$OUT"/m*.txt 2>/dev/null
