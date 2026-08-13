#!/bin/sh
# timing.sh [RESDIR]  -- recover wall-clock cost of a finished run from file mtimes.
#
# run_n.sh / run_n11_P.sh write one output file per bucket, in order, so the gap
# between consecutive mtimes is that bucket's wall time.  Works after the fact:
# nothing needs to have been instrumented.  Multiply by the thread count for a
# CPU-time upper bound (see README "Reporting cost").
RES=${1:-res11}
ls -1 "$RES"/m*.txt 2>/dev/null | head -1 >/dev/null || { echo "no bucket files in $RES"; exit 1; }
if [ "$(uname)" = Darwin ]; then
  mt() { stat -f '%m %N' "$@"; }
else
  mt() { stat -c '%Y %n' "$@"; }
fi
mt "$RES"/m*.txt | sort -n | awk '
  { t=$1; f=$2; sub(/.*\//,"",f)
    if (prev != "") { d=t-prev; printf "%-16s %10.1f s\n", pf, d; tot+=d; if(d>mx){mx=d;mf=pf} }
    prev=t; pf=f }
  END { printf "\n%-16s %10.1f s  = %.2f h  (excludes the first bucket)\n","TOTAL WALL",tot,tot/3600
        printf "%-16s %10.1f s  (%s)\n","largest bucket",mx,mf }'
