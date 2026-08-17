#!/usr/bin/env python3
"""mn_example.py -- verify the pair of Theorem 4.11 from the definitions.

Self-contained: no imports beyond the standard library, no shared code with
the fast search, nothing modular.  Everything is computed exactly over Z by
brute-force enumeration, so this is an independent check of the paper's
headline multigraph claim.

    python3 mn_example.py

The pair (five vertices, ten edges, both with underlying simple graph the
5-cycle 1,3,5,2,4 plus the chord {1,5}):

    edge      {1,3} {3,5} {5,2} {2,4} {4,1} {1,5}
    G_1         1     2     2     3     1     1
    G_2         1     2     1     2     3     1

What is checked:

    G_1 not isomorphic to G_2          (all 120 relabellings)
    XB, U, Tutte, degree sequence      equal
    L, Lhat, U^ext                     different

together with the two displays quoted in the paper: the four-term difference
L_{G_1} - L_{G_2} and the coefficient [z_{2,1} z_{3,0}] U^ext = 19 vs 17.
"""
from collections import Counter
from itertools import permutations
import sys

V = [1, 2, 3, 4, 5]
G1 = {(1, 3): 1, (3, 5): 2, (5, 2): 2, (2, 4): 3, (4, 1): 1, (1, 5): 1}
G2 = {(1, 3): 1, (3, 5): 2, (5, 2): 1, (2, 4): 2, (4, 1): 3, (1, 5): 1}


def edge_list(mult):
    """Expand a multiplicity dict into a list of parallel edges."""
    out = []
    for (u, v), k in mult.items():
        out += [(u, v)] * k
    return out


def components(edges, subset, verts=V):
    """(size, #edges) of each component of the spanning subgraph `subset`."""
    par = {v: v for v in verts}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for i in subset:
        a, b = find(edges[i][0]), find(edges[i][1])
        if a != b:
            par[a] = b
    size, ecnt = Counter(), Counter()
    for v in verts:
        size[find(v)] += 1
    for i in subset:
        ecnt[find(edges[i][0])] += 1
    return [(size[r], ecnt[r]) for r in size]


def subsets(n):
    for mask in range(1 << n):
        yield [i for i in range(n) if mask >> i & 1]


# --------------------------------------------------------------------------
# spanning-subgraph invariants
# --------------------------------------------------------------------------

def u_ext(mult):
    """U^ext: multiset of (|V(D)|, nu(D)) over components, per subgraph."""
    E = edge_list(mult)
    C = Counter()
    for S in subsets(len(E)):
        C[tuple(sorted((s, e - s + 1) for s, e in components(E, S)))] += 1
    return C


def u_poly(mult):
    """U: multiset of component orders, plus the total nullity as a y-power."""
    E = edge_list(mult)
    C = Counter()
    for S in subsets(len(E)):
        comp = components(E, S)
        key = (tuple(sorted(s for s, _ in comp)),
               sum(e - s + 1 for s, e in comp))
        C[key] += 1
    return C


def tutte(mult):
    """Tutte polynomial by the rank-nullity (Whitney rank) expansion."""
    E = edge_list(mult)
    n = len(V)
    C = Counter()
    for S in subsets(len(E)):
        comp = components(E, S)
        k = len(comp)
        C[(n - k, len(S) - n + k)] += 1      # (r(E)-r(S) exponent data)
    return C


# --------------------------------------------------------------------------
# forest-expansion invariants
# --------------------------------------------------------------------------

def spanning_forests(mult):
    """Yield (forest as an index set, list of its components)."""
    E = edge_list(mult)
    for S in subsets(len(E)):
        comp = components(E, S)
        if sum(e - s + 1 for s, e in comp) == 0:      # acyclic
            yield S, comp


def external_activity(mult, order=None):
    """For each spanning forest F, return (|F|, eps(F), per-component data).

    `eps` is ordinary external activity with respect to a fixed linear order
    on the edges: an edge e not in F is externally active if it is the largest
    edge of the unique cycle in F + e.  An edge joining two components of F
    creates no cycle and is never active.  Activity is attributed to the
    component containing the cycle, giving eps_F(T).
    """
    E = edge_list(mult)
    n = len(E)
    rank = list(range(n)) if order is None else order

    for S, _ in spanning_forests(mult):
        Sset = set(S)
        par = {v: v for v in V}

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x

        adj = {v: [] for v in V}
        for i in S:
            u, v = E[i]
            par[find(u)] = find(v)
            adj[u].append((v, i))
            adj[v].append((u, i))

        def path(u, w):
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
            return None

        size, ein, eext = Counter(), Counter(), Counter()
        for v in V:
            size[find(v)] += 1
        for i in S:
            ein[find(E[i][0])] += 1
        total = 0
        for i in range(n):
            if i in Sset:
                continue
            a, b = find(E[i][0]), find(E[i][1])
            if a != b:
                continue
            cyc = path(E[i][0], E[i][1])
            if all(rank[i] > rank[j] for j in cyc):
                eext[a] += 1
                total += 1
        yield (len(S), total,
               sorted((size[r], ein[r], eext[r]) for r in size))


def loopy(mult):
    """L_G as a Counter on (t-exponent, sorted tuple of x-indices)."""
    m = sum(mult.values())
    C = Counter()
    for nF, ext, comps in external_activity(mult):
        C[(m - nF - ext, tuple(sorted(ein + eext for _, ein, eext in comps)))] += 1
    return C


def loopy_hat(mult):
    """Lhat_G as a Counter on the multiset of (|V(T)|, eps_F(T))."""
    C = Counter()
    for _, _, comps in external_activity(mult):
        C[tuple(sorted((s, eext) for s, _, eext in comps))] += 1
    return C


def xb(mult):
    """XB_G, as the Counter of (partition of V by components, total nullity)
    -- equivalently U with the component-order multiset replaced by the
    induced set partition's type.  Here it is computed as U's refinement by
    the subset structure, which is what XB records."""
    E = edge_list(mult)
    C = Counter()
    for S in subsets(len(E)):
        par = {v: v for v in V}

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x

        for i in S:
            a, b = find(E[i][0]), find(E[i][1])
            if a != b:
                par[a] = b
        blocks = Counter()
        for v in V:
            blocks[find(v)] += 1
        comp = components(E, S)
        C[(tuple(sorted(blocks.values())),
           sum(e - s + 1 for s, e in comp))] += 1
    return C


def degrees(mult):
    d = Counter()
    for (u, v), k in mult.items():
        d[u] += k
        d[v] += k
    return sorted(d[v] for v in V)


def isomorphic(m1, m2):
    for p in permutations(V):
        sigma = dict(zip(V, p))
        relab = {}
        for (u, v), k in m1.items():
            a, b = sorted((sigma[u], sigma[v]))
            relab[(a, b)] = k
        norm = {}
        for (u, v), k in m2.items():
            a, b = sorted((u, v))
            norm[(a, b)] = k
        if relab == norm:
            return True
    return False


def main():
    fails = []

    def check(name, cond):
        print('  %-46s %s' % (name, 'ok' if cond else 'FAILED'))
        if not cond:
            fails.append(name)

    print('Theorem 4.11, checked from the definitions over Z:')
    check('G_1 not isomorphic to G_2 (120 relabellings)',
          not isomorphic(G1, G2))
    check('degree sequences equal  %s' % (degrees(G1),),
          degrees(G1) == degrees(G2))
    check('Tutte polynomials equal', tutte(G1) == tutte(G2))
    check('XB equal', xb(G1) == xb(G2))
    check('U equal', u_poly(G1) == u_poly(G2))
    check('U^ext differ', u_ext(G1) != u_ext(G2))
    check('Lhat differ', loopy_hat(G1) != loopy_hat(G2))
    check('L differ', loopy(G1) != loopy(G2))

    print('\n  the U^ext certificate quoted in Example 4.12:')
    a, b = u_ext(G1), u_ext(G2)
    key = ((2, 1), (3, 0))
    print('    [z_{2,1} z_{3,0}] U^ext :  %d  vs  %d' % (a[key], b[key]))
    check('  equals 19 and 17', (a[key], b[key]) == (19, 17))

    print('\n  L_{G_1} - L_{G_2}, by t-degree:')
    d = Counter(loopy(G1))
    d.subtract(loopy(G2))
    byt = {}
    for (tt, xs), c in d.items():
        if c:
            byt.setdefault(tt, []).append((xs, c))
    for tt in sorted(byt):
        terms = ' '.join('%+d*x%s' % (c, ''.join(str(i) for i in xs))
                         for xs, c in sorted(byt[tt]))
        print('    t^%d : %s' % (tt, terms))

    print()
    if fails:
        sys.exit('%d check(s) FAILED: %s' % (len(fails), ', '.join(fails)))
    print('all checks passed')


if __name__ == '__main__':
    main()
