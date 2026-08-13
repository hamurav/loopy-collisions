#!/usr/bin/env python3
"""
verdict.py -- aggregate the CLASS lines of one or more loopy_collisions runs
(e.g. the per-edge-count jobs of run_n.sh) and decide whether the three
partitions coincide.

    python3 verdict.py res10/*.txt

Logic (see README):  Lhat always refines L, so every L-class is a union of
Lhat-classes.  For n <= 10, Markstrom's Observation 6.6 gives Lhat <=> XB.
Hence at n <= 10 the entire content of Conjecture 4.2 is the single line

        L-classes == Lhat-classes ?

and the comparison Lhat vs XB is an independent re-derivation of Observation 6.6.
"""
import sys
import collections


def read_classes(paths):
    cls = {k: set() for k in ("L", "Lhat", "XB")}
    edges = {k: [] for k in cls}
    for p in paths:
        for line in open(p):
            if not line.startswith("CLASS"):
                continue
            parts = line.split()
            inv = parts[1]
            m = int(parts[3].split('=')[1])
            members = frozenset(parts[5:])
            if inv in cls:
                cls[inv].add(members)
                edges[inv].append(m)
    return cls, edges


def main():
    paths = sys.argv[1:]
    if not paths:
        print("usage: verdict.py <loopy_collisions output files>")
        return 1
    cls, edges = read_classes(paths)

    print("class counts")
    for k in ("L", "Lhat", "XB"):
        e = edges[k]
        rng = f"{min(e)}..{max(e)}" if e else "-"
        sizes = collections.Counter(len(c) for c in cls[k])
        print(f"  {k:5s} {len(cls[k]):6d} classes   edges {rng:10s} sizes {dict(sizes)}")

    # Lhat must refine L; check it explicitly rather than assuming
    refines = all(any(c <= C for C in cls["L"]) for c in cls["Lhat"])
    print(f"\nsanity: every Lhat-class sits inside an L-class : "
          f"{'yes' if refines else 'NO  (something is wrong)'}")

    LH = cls["L"] == cls["Lhat"]
    HX = cls["Lhat"] == cls["XB"]
    LX = cls["L"] == cls["XB"]

    print()
    print(f"VERDICT  L == Lhat  : {'YES' if LH else 'NO'}"
          "    <- Conjecture 4.2 at this n (given Markstrom Obs. 6.6)")
    print(f"VERDICT  Lhat == XB : {'YES' if HX else 'NO'}"
          "    <- independent re-derivation of Markstrom Obs. 6.6")
    print(f"VERDICT  L == XB    : {'YES' if LX else 'NO'}")

    if not LH:
        extra = cls["L"] - cls["Lhat"]
        print("\n*** L-classes that are NOT single Lhat-classes "
              "(candidate counterexamples to Conjecture 4.2):")
        for c in sorted(extra, key=sorted):
            print("   ", " ".join(sorted(c)))
    if not HX:
        extra = cls["XB"] - cls["Lhat"]
        print("\n*** XB-classes that split under Lhat "
              "(would answer the Merino-Noble question):")
        for c in sorted(extra, key=sorted):
            print("   ", " ".join(sorted(c)))

    print("\nNow confirm exactly:  cat <files> | python3 exact_check.py")
    return 0 if (LH and HX and LX) else 3


if __name__ == "__main__":
    sys.exit(main())
