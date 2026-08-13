#!/usr/bin/env python3
"""
exact_check.py -- exact (non-randomized) recomputation of the candidate
collision classes produced by loopy_collisions.

The fingerprints in loopy_collisions are evaluations, so fingerprint equality is
a NECESSARY condition for equality of the invariants: no genuine collision can
be missed.  False positives, however, are possible in principle, so every
reported class must be confirmed exactly.  That is what this script does, with
exact integer arithmetic, for

    L(1,x)  : the loopy polynomial            (monomials in x_0,x_1,...)
    Lhat    : the refined loopy polynomial    ( = U^ext ; monomials in u_{s,r} )
    XB      : the Tutte symmetric function    (coefficients in q, basis p_lambda)

Usage
-----
    ./loopy_collisions < c9.g6 | ./exact_check.py
    ./exact_check.py res9.txt
    echo "GCRdfW GCRRVW" | ./exact_check.py --raw

Reads CLASS lines of the form

    CLASS <inv> size=<k> m=<m> : <g6> <g6> ...

(or, with --raw, whitespace-separated groups of graph6 strings, one group per
line) and prints, for each distinct class, whether its members really do share
each of the three invariants.

Cost is O(3^n) polynomial operations per graph, so n = 10, 11 are fine.
"""
import sys
import collections
from itertools import combinations


# ------------------------------------------------------------------ graph6 --
def g6_decode(s):
    n = ord(s[0]) - 63
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        bits += [(v >> k) & 1 for k in range(5, -1, -1)]
    E, p = [], 0
    for j in range(1, n):
        for i in range(j):
            if bits[p]:
                E.append((i, j))
            p += 1
    return n, E


# --------------------------------------------------- connected-block polys --
def blocks(n, E):
    """ec[S] = #edges inside S ; C[S] = list of coefficients of C_G(S;q)."""
    ec = [0] * (1 << n)
    for (i, j) in E:
        bi, bj = 1 << i, 1 << j
        for S in range(1 << n):
            if (S & bi) and (S & bj):
                ec[S] += 1
    maxe = ec[(1 << n) - 1]
    binom = [[0] * (maxe + 2) for _ in range(maxe + 2)]
    for i in range(maxe + 2):
        binom[i][0] = 1
        for j in range(1, i + 1):
            binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j]

    def polmul(a, b):
        c = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        c[i + j] += x * y
        return c

    N = [[binom[ec[S]][k] for k in range(ec[S] + 1)] for S in range(1 << n)]
    C = [None] * (1 << n)
    C[0] = [0]
    for S in range(1, 1 << n):
        low = S & (-S)
        rest = S ^ low
        acc = list(N[S])
        T = rest
        while True:
            T = (T - 1) & rest if T else -1
            if T < 0:
                break
            B = T | low
            p = polmul(C[B], N[S ^ B])
            if len(p) > len(acc):
                acc += [0] * (len(p) - len(acc))
            for i, v in enumerate(p):
                acc[i] -= v
            if T == 0:
                break
        while len(acc) > 1 and acc[-1] == 0:
            acc.pop()
        C[S] = acc
    return ec, C


def to_activity_basis(poly, s):
    """C(q) = q^{s-1} * sum_r a_r (1+q)^r  ->  {r: a_r}.  Returns {} if C = 0."""
    if all(v == 0 for v in poly):
        return {}
    assert all(v == 0 for v in poly[:s - 1]), (poly, s)
    p = poly[s - 1:]
    binom = [[0] * (len(p) + 1) for _ in range(len(p) + 1)]
    for i in range(len(p) + 1):
        binom[i][0] = 1
        for j in range(1, i + 1):
            binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j]
    out = collections.defaultdict(int)
    for i, c in enumerate(p):
        if c:
            for r in range(i + 1):
                out[r] += c * binom[i][r] * ((-1) ** (i - r))
    return {k: v for k, v in out.items() if v}


# ---------------------------------------------------- subset DP over blocks --
def _partition_dp(n, weight):
    """sum over set partitions of prod weight(B), with monomial bookkeeping.

    weight(B) must return a dict {monomial-part: coefficient}; monomial parts
    are combined by tuple concatenation and sorting."""
    Z = [None] * (1 << n)
    Z[0] = {(): 1}
    for S in range(1, 1 << n):
        low = S & (-S)
        rest = S ^ low
        acc = collections.defaultdict(int)
        T = rest
        while True:
            B = T | low
            w = weight(B)
            if w:
                for mono, c in Z[S ^ B].items():
                    for part, d in w.items():
                        acc[tuple(sorted(mono + (part,)))] += c * d
            if T == 0:
                break
            T = (T - 1) & rest
        Z[S] = {k: v for k, v in acc.items() if v}
    return Z[(1 << n) - 1]


def L_exact(n, E):
    """L_G(1,x): dict {sorted tuple of x-indices: coefficient}."""
    ec, C = blocks(n, E)
    cache = {}

    def w(B):
        if B not in cache:
            s = bin(B).count('1')
            cache[B] = {s - 1 + r: a for r, a in to_activity_basis(C[B], s).items()}
        return cache[B]
    return _partition_dp(n, w)


def Lhat_exact(n, E):
    """refined loopy polynomial: dict {sorted tuple of (s,r): coefficient}."""
    ec, C = blocks(n, E)
    cache = {}

    def w(B):
        if B not in cache:
            s = bin(B).count('1')
            cache[B] = {(s, r): a for r, a in to_activity_basis(C[B], s).items()}
        return cache[B]
    return _partition_dp(n, w)


def XB_exact(n, E):
    """XB_G: dict {(sorted tuple of block sizes, exponent of q): coefficient}."""
    ec, C = blocks(n, E)
    Z = [None] * (1 << n)
    Z[0] = {((), 0): 1}
    for S in range(1, 1 << n):
        low = S & (-S)
        rest = S ^ low
        acc = collections.defaultdict(int)
        T = rest
        while True:
            B = T | low
            b = bin(B).count('1')
            if any(C[B]):
                for (mono, deg), c in Z[S ^ B].items():
                    for e, d in enumerate(C[B]):
                        if d:
                            acc[(tuple(sorted(mono + (b,))), deg + e)] += c * d
            if T == 0:
                break
            T = (T - 1) & rest
        Z[S] = {k: v for k, v in acc.items() if v}
    return Z[(1 << n) - 1]


# ------------------------------------------------------------------- driver --
def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    raw = '--raw' in sys.argv
    src = open(args[0]) if args else sys.stdin

    groups = set()
    for line in src:
        line = line.strip()
        if not line:
            continue
        if raw:
            parts = line.split()
        elif line.startswith('CLASS'):
            parts = line.split(':')[1].split() if ':' in line else line.split()[4:]
        else:
            continue
        if len(parts) > 1:
            groups.add(frozenset(parts))

    if not groups:
        print("no candidate classes found on input")
        return

    print(f"{len(groups)} distinct candidate classes to confirm\n")
    stats = collections.Counter()
    for g in sorted(groups, key=lambda s: sorted(s)):
        codes = sorted(g)
        decoded = [g6_decode(c) for c in codes]
        n = decoded[0][0]
        Ls = [L_exact(*d) for d in decoded]
        Hs = [Lhat_exact(*d) for d in decoded]
        Xs = [XB_exact(*d) for d in decoded]
        sameL = all(x == Ls[0] for x in Ls)
        sameH = all(x == Hs[0] for x in Hs)
        sameX = all(x == Xs[0] for x in Xs)
        stats['L'] += sameL
        stats['Lhat'] += sameH
        stats['XB'] += sameX
        stats['total'] += 1
        flag = '' if (sameL and sameH and sameX) else '   <-- NOT all equal'
        print(f"  n={n} m={len(decoded[0][1])}  {' '.join(codes)}   "
              f"L:{'=' if sameL else '#'} Lhat:{'=' if sameH else '#'} "
              f"XB:{'=' if sameX else '#'}{flag}")

    t = stats['total']
    print(f"\nconfirmed exactly:  L {stats['L']}/{t}   Lhat {stats['Lhat']}/{t}   "
          f"XB {stats['XB']}/{t}")
    if stats['L'] == stats['Lhat'] == stats['XB'] == t:
        print("all candidate classes are genuine, and identical for the three invariants")


if __name__ == '__main__':
    main()
