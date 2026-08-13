#!/usr/bin/env python3
"""make_lists.py RESDIR [-o OUT] [--per-invariant]

Consolidate the per-bucket output of a collision run into one canonical list.

A run leaves its results scattered over one file per edge bucket (and per shard
of the sharded buckets), each containing CLASS lines for all three invariants.
This collects them, removes duplicates, checks that the three invariants agree,
sorts everything deterministically, and writes a single documented file.

OUTPUT FORMAT  one class per line, blank- and comment-lines ignored by the
other tools in this directory:

    n  m  g6 g6 [g6 ...]

n = vertices, m = edges, then the members in graph6, sorted.  The class size is
the number of graph6 fields, so classes of size 2, 3, 4, ... all fit the same
format and are read correctly by comp_check.c and lift.py.

    python3 make_lists.py res11 -o data/collisions_n11.txt
    python3 make_lists.py res10 -o data/collisions_n10.txt

With --per-invariant, additionally writes OUT.L, OUT.Lhat and OUT.XB.  That is
only useful if the invariants ever disagree; if they do, this script says so
loudly and exits nonzero, because the whole point of the run is that they do
not.
"""
import sys, os, glob, argparse
from collections import defaultdict, Counter

INVARIANTS = ("L", "Lhat", "XB")


def read_classes(resdir):
    """-> {invariant: {(m, members) }},  n,  number of files read"""
    files = sorted(glob.glob(os.path.join(resdir, "*.txt"))) \
        if os.path.isdir(resdir) else [resdir]
    if not files:
        sys.exit("make_lists: no .txt files in %s" % resdir)
    out = {inv: set() for inv in INVARIANTS}
    n = None
    for f in files:
        for line in open(f):
            if line.startswith("SUMMARY"):
                nn = int(line.split()[2].split("=")[1])
                if n is None:
                    n = nn
                elif nn != n:
                    sys.exit("make_lists: mixed vertex counts (%d and %d) in %s"
                             % (n, nn, f))
            elif line.startswith("CLASS"):
                p = line.split()
                inv = p[1]
                m = int(p[3].split("=")[1])
                members = tuple(sorted(line.split(":", 1)[1].split()))
                if inv in out:
                    out[inv].add((m, members))
    if n is None:
        sys.exit("make_lists: no SUMMARY line found -- is %s a results dir?" % resdir)
    return out, n, len(files)


def write(path, n, classes, note):
    classes = sorted(classes, key=lambda c: (c[0], len(c[1]), c[1]))
    sizes = Counter(len(mem) for _, mem in classes)
    ms = [m for m, _ in classes]
    with open(path, "w") as fh:
        fh.write("# Collision classes of the loopy polynomial L_G.\n")
        fh.write("# %s\n" % note)
        fh.write("#\n")
        fh.write("# Non-isomorphic connected simple graphs on n = %d vertices with equal\n" % n)
        fh.write("# L_G.  In this range the classes of the refined loopy polynomial Lhat_G\n")
        fh.write("# and of the Tutte symmetric function XB_G are identical.\n")
        fh.write("#\n")
        fh.write("# Format:  n  m  graph6 [graph6 ...]   (class size = number of graph6 fields)\n")
        fh.write("#\n")
        if classes:
            fh.write("# classes    : %d\n" % len(classes))
            fh.write("# by size    : %s\n" % ", ".join(
                "%d of size %d" % (k, s) for s, k in sorted(sizes.items())))
            fh.write("# graphs     : %d\n" % sum(s * k for s, k in sizes.items()))
            fh.write("# pairs      : %d\n" % sum(s * (s - 1) // 2 * k
                                                 for s, k in sizes.items()))
            fh.write("# edge range : %d..%d\n" % (min(ms), max(ms)))
        else:
            fh.write("# classes    : 0\n")
        fh.write("#\n")
        for m, mem in classes:
            fh.write("%d %d %s\n" % (n, m, " ".join(mem)))
    return len(classes), sizes, ms


def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("resdir")
    ap.add_argument("-o", "--out")
    ap.add_argument("--per-invariant", action="store_true")
    ap.add_argument("-h", "--help", action="store_true")
    a = ap.parse_args()
    if a.help:
        print(__doc__)
        return 0

    per_inv, n, nfiles = read_classes(a.resdir)
    out = a.out or "collisions_n%d.txt" % n

    # ---- agreement -------------------------------------------------------
    base = per_inv["L"]
    disagree = False
    for inv in INVARIANTS[1:]:
        if not per_inv[inv]:
            print("warning: no %s classes found" % inv)
            continue
        only_l, only_i = base - per_inv[inv], per_inv[inv] - base
        if only_l or only_i:
            disagree = True
            print("!! L and %s DISAGREE: %d only in L, %d only in %s"
                  % (inv, len(only_l), len(only_i), inv))
            for m, mem in sorted(only_l)[:5]:
                print("     only L    m=%d  %s" % (m, " ".join(mem)))
            for m, mem in sorted(only_i)[:5]:
                print("     only %-4s m=%d  %s" % (inv, m, " ".join(mem)))

    note = ("Read from %d result file(s) in %s.  The L-, Lhat- and XB-classes "
            "coincide." % (nfiles, a.resdir)) if not disagree else \
           ("Read from %d result file(s) in %s.  WARNING: the invariants do "
            "NOT all agree." % (nfiles, a.resdir))

    ncls, sizes, ms = write(out, n, base, note)
    print("wrote %s" % out)
    print("  n = %d, %d classes from %d file(s)" % (n, ncls, nfiles))
    if ncls:
        print("  sizes: %s" % ", ".join("%d x size %d" % (k, s)
                                        for s, k in sorted(sizes.items())))
        print("  graphs: %d, unordered pairs: %d, edges %d..%d"
              % (sum(s * k for s, k in sizes.items()),
                 sum(s * (s - 1) // 2 * k for s, k in sizes.items()),
                 min(ms), max(ms)))

    if a.per_invariant:
        for inv in INVARIANTS:
            if per_inv[inv]:
                write("%s.%s" % (out, inv), n, per_inv[inv],
                      "Classes of %s only." % inv)
                print("  wrote %s.%s" % (out, inv))

    if disagree:
        print("\nthe invariants disagree -- do not publish this list as it stands")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
