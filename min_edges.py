#!/usr/bin/env python3
"""min_edges.py -- how few edges can a multigraph counterexample have?

Exhaustive search over connected loopless multigraphs with a bounded number of
EDGES (not, as in the main search, a bounded multiplicity).  It looks for the
pairs that matter for Theorem 4.13 of the paper:

    non-isomorphic G, H  with  XB_G = XB_H  but  L_G != L_H.

    ./gen_graphs 8 | python3 min_edges.py --maxm 9
    geng -q -c 8      | python3 min_edges.py --maxm 9        # with nauty

Input is connected simple graphs in graph6 on stdin; these are read as the
possible *underlying* simple graphs.  For each one with at most `maxm` edges,
every assignment of positive multiplicities with total at most `maxm` is
enumerated, one representative per orbit of Aut(G) acting on E(G), so the
multigraphs tested are pairwise non-isomorphic.  L and XB are then computed
exactly over Z from their defining expansions -- no fingerprints, no modular
arithmetic, nothing shared with the fast C code.

WHY THIS IS A COMPLETE ANSWER FOR m <= 9

A connected graph with m edges has at most m+1 vertices, so a counterexample
with at most eight edges would live on at most nine vertices.  Running this
script for n = 2,...,9 therefore settles the question for m <= 8 outright.
For m = 9 the same runs settle n <= 9, and n = 10 is immediate: a connected
graph on ten vertices with nine edges is a tree, hence simple, and forests are
separated by all of these invariants (Corollary 3.19 of the paper).

RESULT (the runs reported in Section 4.4.2)

    n     multigraphs with m <= 9     collisions
    2                    9            none
    3                   43            none
    4                  233            none
    5                  722            none
    6                 1462            one, at m = 9
    7                 1738            one, at m = 9
    8                 1275            one, at m = 9
    9                  526            one, at m = 9

So no counterexample has fewer than nine edges, and nine is attained.  The
six-vertex witness has maximum multiplicity two and therefore also appears in
the (6,2) row of the table in Section 4.4.2; it is the single XB-class in the
m = 9 bucket of the (6,4) run.

Cost is dominated by computing Aut(G) by brute force over the n! vertex
permutations, so n = 9 takes a couple of minutes and n = 10 is not attempted
(it is not needed -- see above).
"""
import argparse
import sys
from collections import Counter
from itertools import permutations


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------

def read_graph6(line):
    """decode one graph6 string into (n, sorted list of edges)"""
    s = line.strip()
    n = ord(s[0]) - 63
    bits = []
    for c in s[1:]:
        v = ord(c) - 63
        bits += [(v >> k) & 1 for k in range(5, -1, -1)]
    edges, p = [], 0
    for j in range(1, n):
        for i in range(j):
            if bits[p]:
                edges.append((i, j))
            p += 1
    return n, edges


def automorphisms(n, edges):
    """Aut(G) as a list of vertex permutations, by brute force."""
    E = set(edges)
    out = []
    for pm in permutations(range(n)):
        if {tuple(sorted((pm[u], pm[v]))) for u, v in edges} == E:
            out.append(pm)
    return out


def multiplicity_orbits(n, edges, maxm):
    """One multiplicity vector per Aut(G)-orbit, with total at most maxm.

    Multiplicities are at least one, so the underlying simple graph is
    exactly `edges` and the resulting multigraphs are pairwise
    non-isomorphic."""
    k = len(edges)
    if k == 0 or k > maxm:
        return []
    index = {e: i for i, e in enumerate(edges)}
    maps = [[index[tuple(sorted((pm[u], pm[v])))] for (u, v) in edges]
            for pm in automorphisms(n, edges)]
    out = []

    def rec(i, remaining, cur):
        if i == k:
            t = tuple(cur)
            if min(tuple(t[m[j]] for j in range(k)) for m in maps) == t:
                out.append(t)
            return
        # leave at least one unit for each remaining edge
        for c in range(1, remaining - (k - i - 1) + 1):
            cur.append(c)
            rec(i + 1, remaining - c, cur)
            cur.pop()

    rec(0, maxm, [])
    return out


# --------------------------------------------------------------------------
# the two invariants, computed exactly from their definitions
# --------------------------------------------------------------------------

def components(n, E, subset):
    """(order, edge count) of each component of the spanning subgraph."""
    par = list(range(n))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for i in subset:
        a, b = find(E[i][0]), find(E[i][1])
        if a != b:
            par[a] = b
    order, ecount = Counter(), Counter()
    for v in range(n):
        order[find(v)] += 1
    for i in subset:
        ecount[find(E[i][0])] += 1
    return [(order[r], ecount[r]) for r in order]


def XB(n, E):
    """XB_G as the multiset of (|A|, multiset of component orders) over all
    spanning subgraphs A -- the expansion XB_G = sum_A q^|A| prod_D p_{|V(D)|}."""
    C = Counter()
    for mask in range(1 << len(E)):
        S = [i for i in range(len(E)) if mask >> i & 1]
        C[(len(S), tuple(sorted(s for s, _ in components(n, E, S))))] += 1
    return frozenset(C.items())


def L(n, E):
    """L_G as the multiset of (t-exponent, multiset of x-indices), from the
    spanning-forest expansion with external activity in the edge order given
    by the position in E."""
    m = len(E)
    C = Counter()
    for mask in range(1 << m):
        S = [i for i in range(m) if mask >> i & 1]
        comp = components(n, E, S)
        if sum(e - s + 1 for s, e in comp):          # not acyclic
            continue
        par = list(range(n))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x

        adj = {v: [] for v in range(n)}
        for i in S:
            u, v = E[i]
            par[find(u)] = find(v)
            adj[u].append((v, i))
            adj[v].append((u, i))

        def tree_path(u, w):
            """edge indices on the F-path from u to w"""
            stack, seen = [(u, [])], {u}
            while stack:
                x, p = stack.pop()
                if x == w:
                    return p
                for y, i in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append((y, p + [i]))
            return []

        order, inside, active = Counter(), Counter(), Counter()
        for v in range(n):
            order[find(v)] += 1
        for i in S:
            inside[find(E[i][0])] += 1
        chosen, total = set(S), 0
        for i in range(m):
            if i in chosen:
                continue
            a, b = find(E[i][0]), find(E[i][1])
            if a != b:                                # joins two components
                continue
            if all(i > j for j in tree_path(E[i][0], E[i][1])):
                active[a] += 1
                total += 1
        C[(m - len(S) - total,
           tuple(sorted(inside[r] + active[r] for r in order)))] += 1
    return frozenset(C.items())


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--maxm', type=int, default=9,
                    help='maximum total number of edges (default 9)')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='print every multigraph tested')
    args = ap.parse_args()

    total, by_class, n_seen = 0, {}, set()
    for line in sys.stdin:
        if not line.strip() or line.startswith('>'):
            continue
        n, edges = read_graph6(line)
        n_seen.add(n)
        if len(edges) > args.maxm:
            continue
        for mv in multiplicity_orbits(n, edges, args.maxm):
            E = []
            for c, e in zip(mv, edges):
                E += [e] * c
            total += 1
            if args.verbose:
                print('  n=%d m=%d %s %s' % (n, len(E), edges, mv))
            by_class.setdefault((n, len(E), XB(n, E)), []).append((edges, mv, L(n, E)))

    bad = [(k[0], k[1], v) for k, v in by_class.items()
           if len(v) > 1 and len({x[2] for x in v}) > 1]

    label = ','.join(str(x) for x in sorted(n_seen))
    print('n = %s, m <= %d: %d connected loopless multigraphs tested'
          % (label, args.maxm, total))
    if not bad:
        print('  no pair with equal XB and different L')
        return
    print('  %d XB-class(es) split by L:' % len(bad))
    for n, m, v in sorted(bad):
        print('    n=%d m=%d' % (n, m))
        for edges, mv, _ in v:
            print('      underlying %s  multiplicities %s'
                  % (' '.join('%d-%d' % e for e in edges),
                     ' '.join(map(str, mv))))


if __name__ == '__main__':
    main()
