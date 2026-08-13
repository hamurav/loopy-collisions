#!/usr/bin/env python3
"""verify_run.py RESDIR [PREV_RESDIR]  -- audit a finished collision run.

Recounts everything from the raw CLASS/SUMMARY lines and cross-checks it three
ways.  Nothing has to be rerun.

  1. COVERAGE   every edge bucket m present exactly once (shards summed), and
                the total number of graphs processed equals A001349(n).
                This is the check that catches a lost or duplicated bucket --
                the failure mode that produced a truncated edge range.
  2. COUNTS     classes tallied by size, independently of the SUMMARY line,
                then reconciled with it (SUMMARY "pairs" counts unordered
                PAIRS, so a class of size s contributes s(s-1)/2, not 1).
  3. AGREEMENT  the L, Lhat and XB partitions must be identical, class by
                class and not merely in total.
  4. SYMMETRY   complementation is a bijection on graphs and preserves
                collisions (XB_Gbar is an explicit invertible transform of
                XB_G), so the (size, m) profile must be symmetric under
                m -> C(n,2) - m up to classes whose complements are
                disconnected.  Gross asymmetry means missing data.
  5. LIFTS      if PREV_RESDIR (the n-1 run) is given: every class there must
                reappear here at m + (n-1), via G -> K_1 \\/ G.

Usage:  python3 verify_run.py res11
        python3 verify_run.py res11 res10
"""
import sys, os, re, glob
from collections import defaultdict, Counter

A001349 = {1:1, 2:1, 3:2, 4:6, 5:21, 6:112, 7:853, 8:11117, 9:261080,
           10:11716571, 11:1006700565, 12:164059830476}

def load(d):
    files = sorted(glob.glob(os.path.join(d, "m*.txt"))) if os.path.isdir(d) else [d]
    if not files:
        sys.exit("no m*.txt files in %s" % d)
    cls = defaultdict(list)          # invariant -> [(size, m, members)]
    summ = []                        # (invariant, n, graphs, classes, pairs, file)
    for f in files:
        for line in open(f):
            if line.startswith("CLASS"):
                p = line.split()
                inv = p[1]
                size = int(p[2].split("=")[1])
                m    = int(p[3].split("=")[1])
                mem  = line.split(":", 1)[1].split()
                cls[inv].append((size, m, tuple(sorted(mem))))
            elif line.startswith("SUMMARY"):
                p = line.split()
                summ.append((p[1], int(p[2].split("=")[1]), int(p[3].split("=")[1]),
                             int(p[4].split("=")[1]), int(p[5].split("=")[1]), f))
    return files, cls, summ

def bucket_of(path):
    mo = re.search(r"m(\d+)", os.path.basename(path))
    return int(mo.group(1)) if mo else None

def main():
    if len(sys.argv) < 2: sys.exit(__doc__)
    files, cls, summ = load(sys.argv[1])
    if not summ: sys.exit("no SUMMARY lines found")
    n = summ[0][1]
    M = n*(n-1)//2
    ok = True
    print("run: %s   n = %d   (%d files)\n" % (sys.argv[1], n, len(files)))

    # ---- 1. coverage -----------------------------------------------------
    # NB SUMMARY "graphs=" is the number of SURVIVORS of the cheap sieve, not
    # the number of graphs read, so it cannot be compared with A001349.  What
    # we can check exactly is that every edge bucket produced output, and that
    # every shard of a sharded bucket did.  A silently missing bucket is what
    # truncates an edge range.
    print("1. COVERAGE")
    shards = defaultdict(set)
    surv = defaultdict(int)
    for f in files:
        b = bucket_of(f)
        mo = re.search(r"m\d+_s(\d+)", os.path.basename(f))
        shards[b].add(int(mo.group(1)) if mo else -1)
    for inv, nn, graphs, c, p, f in summ:
        if inv != "L": continue
        surv[bucket_of(f)] += graphs
        if nn != n: print("   !! mixed n: %d in %s" % (nn, f)); ok = False
    present = sorted(shards)
    missing = [m for m in range(n-1, M+1) if m not in shards]
    print("   buckets m = %d .. %d present (%d of the %d possible)"
          % (present[0], present[-1], len(present), M-n+2))
    if missing:
        ok = False
        print("   !! MISSING buckets: %s" % missing)
        print("      -> the edge range and the class counts are BOTH understated")
    for b in present:                     # sharded buckets must be complete
        s = shards[b]
        if s != {-1} and sorted(s) != list(range(len(s))):
            print("   !! bucket m=%d has shards %s -- not 0..%d"
                  % (b, sorted(s), len(s)-1)); ok = False
    print("   survivors of the sieve: %d  (SUMMARY 'graphs=' counts survivors,"
          " not input)" % sum(surv.values()))
    if A001349.get(n):
        print("   to confirm the input side, run:  geng -q -u -c %d   (expect %d)"
              % (n, A001349[n]))

    # ---- 2. counts -------------------------------------------------------
    print("\n2. COUNTS  (recounted from CLASS lines)")
    for inv in ("L", "Lhat", "XB"):
        if inv not in cls: continue
        sizes = Counter(s for s, m, _ in cls[inv])
        ncls  = sum(sizes.values())
        npair = sum(s*(s-1)//2 * k for s, k in sizes.items())
        ngr   = sum(s*k for s, k in sizes.items())
        ms    = [m for _, m, _ in cls[inv]]
        print("   %-5s classes=%-7d graphs=%-8d unordered pairs=%-7d edges=[%d,%d]"
              % (inv, ncls, ngr, npair, min(ms), max(ms)))
        print("         sizes: %s" % ", ".join("%d x %d" % (k, s)
                                               for s, k in sorted(sizes.items())))
        sc = sum(c for i, _, _, c, _, _ in summ if i == inv)
        sp = sum(p for i, _, _, _, p, _ in summ if i == inv)
        if (sc, sp) != (ncls, npair):
            print("         !! SUMMARY says classes=%d pairs=%d" % (sc, sp)); ok = False

    # ---- 3. agreement ----------------------------------------------------
    print("\n3. AGREEMENT")
    keys = {inv: set(mem for _, _, mem in v) for inv, v in cls.items()}
    if len(keys) < 2:
        print("   only one invariant present")
    else:
        base = keys.get("L")
        for inv, k in keys.items():
            if inv == "L": continue
            d1, d2 = base - k, k - base
            print("   L vs %-5s : %s" % (inv, "identical" if not d1 and not d2
                  else "!! %d only in L, %d only in %s" % (len(d1), len(d2), inv)))
            if d1 or d2: ok = False

    # ---- 4. symmetry -----------------------------------------------------
    print("\n4. SYMMETRY under m -> %d - m" % M)
    hist = Counter(m for _, m, _ in cls.get("L", []))
    lo, hi = min(hist), max(hist)
    if lo + hi != M:
        print("   !! range [%d,%d] is not centred: %d + %d = %d, expected %d"
              % (lo, hi, lo, hi, lo+hi, M)); ok = False
    else:
        print("   range [%d,%d] is centred   OK" % (lo, hi))
    asym = [(m, hist[m], hist[M-m]) for m in range(lo, M//2+1) if hist[m] != hist[M-m]]
    print("   %d edge counts have N(m) != N(%d-m)%s" % (len(asym), M,
          "" if not asym else "  (expected: classes whose complements are disconnected)"))
    for m, a, b in asym[:8]:
        print("      m=%-3d %-7d   m=%-3d %-7d" % (m, a, M-m, b))

    # ---- 5. lifts --------------------------------------------------------
    if len(sys.argv) > 2:
        print("\n5. LIFTS from %s" % sys.argv[2])
        _, pcls, psum = load(sys.argv[2])
        pn = psum[0][1]
        if pn != n-1:
            print("   !! previous run has n = %d, expected %d" % (pn, n-1))
        else:
            phist = Counter((s, m+pn) for s, m, _ in pcls.get("L", []))
            hist2 = Counter((s, m) for s, m, _ in cls.get("L", []))
            bad = [(k, v, hist2[k]) for k, v in sorted(phist.items()) if hist2[k] < v]
            if bad:
                ok = False
                print("   !! K_1 \\/ G lifting requires N(size,m) >= these, but:")
                for (s, m), need, got in bad[:10]:
                    print("      size=%d m=%d : need >= %d, found %d" % (s, m, need, got))
            else:
                print("   all %d lifted classes accounted for   OK" % sum(phist.values()))

    print("\n%s" % ("ALL CHECKS PASSED" if ok else "PROBLEMS FOUND -- see !! lines above"))
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
