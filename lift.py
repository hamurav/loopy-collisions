#!/usr/bin/env python3
"""lift.py < classes.txt > lifted.txt

Apply G -> K_1 \\/ G (join a new vertex to every existing one) to each class.

Since the complement of K_1 \\/ G is K_1 (disjoint union) complement of G, and
XB_Gbar is an explicit invertible transform of XB_G, a collision class on n
vertices with m edges maps to a collision class on n+1 vertices with m+n edges,
of the same size.  Piping a run's list through this and grepping the next run's
list is a check on the search that does not go through the fingerprints at all.

Reads and writes the make_lists.py format:  n  m  g6 [g6 ...]
Classes of any size are handled.
"""
import sys


def decode(s):
    n = ord(s[0]) - 63
    adj = [[0] * n for _ in range(n)]
    p, c = 0, s[1:]
    for j in range(1, n):
        for i in range(j):
            if (ord(c[p // 6]) - 63) >> (5 - (p % 6)) & 1:
                adj[i][j] = adj[j][i] = 1
            p += 1
    return n, adj


def encode(n, adj):
    bits = "".join(str(adj[i][j]) for j in range(1, n) for i in range(j))
    bits += "0" * (-len(bits) % 6)
    return chr(n + 63) + "".join(chr(int(bits[k:k + 6], 2) + 63)
                                for k in range(0, len(bits), 6))


def join_one(n, adj):
    """K_1 \\/ G : new vertex n adjacent to all of 0..n-1"""
    m = n + 1
    a = [[0] * m for _ in range(m)]
    for i in range(n):
        for j in range(n):
            a[i][j] = adj[i][j]
        a[i][n] = a[n][i] = 1
    return m, a


def main():
    for line in sys.stdin:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        n, m, members = int(f[0]), int(f[1]), f[2:]
        if len(members) < 2:
            continue
        lifted = [encode(*join_one(*decode(g))) for g in members]
        print(n + 1, m + n, " ".join(sorted(lifted)))


if __name__ == "__main__":
    main()
