/*  gen_multi.c  --  generate loopless multigraphs, in multig -T text format
 *
 *  A self-contained stand-in for nauty's `multig`, so that mg_collisions.c can
 *  be exercised without nauty.  If you have nauty, prefer
 *
 *      geng -q -c N | multig -mMU -T | ./mg_collisions
 *
 *  which is faster and better tested.  This program does the same job:
 *
 *      ./gen_graphs N | ./gen_multi MU | ./mg_collisions
 *
 *  METHOD
 *    Read connected simple graphs in graph6, one per line.  Two multigraphs
 *    with the same underlying simple graph G are isomorphic exactly when their
 *    multiplicity vectors lie in the same orbit of Aut(G) acting on E(G).  So
 *    for each G we compute Aut(G) by brute force over the n! vertex
 *    permutations, then enumerate multiplicity vectors in {1,...,MU}^{E(G)} and
 *    keep the lexicographically least member of each orbit.  The output is
 *    therefore pairwise non-isomorphic, as multig's is.
 *
 *    Cost is MU^{|E|} per input graph, so the dense graphs dominate: at n = 7
 *    the complete graph alone contributes MU^{21}.
 *
 *  BUILD  cc -O3 -o gen_multi gen_multi.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 9
#define MAXE (MAXN*(MAXN-1)/2)
#define MAXAUT 362880

static int n, nedge;
static int eu[MAXE], ev[MAXE];          /* edges of the underlying graph      */
static int pidx[MAXN][MAXN];            /* vertex pair -> edge index or -1    */
static int perm[MAXN], used[MAXN];
static int aut[MAXAUT][MAXE];           /* Aut(G) as permutations of E(G)     */
static int naut;
static int adj[MAXN];

static int g6_decode(const char *s, int *a)
{
    int nn = (int)s[0] - 63;
    if (nn < 1 || nn > MAXN) return -1;
    for (int i = 0; i < nn; i++) a[i] = 0;
    int p = 0; const char *c = s + 1;
    for (int j = 1; j < nn; j++)
        for (int i = 0; i < j; i++, p++)
            if (((int)c[p/6] - 63) >> (5 - (p%6)) & 1) { a[i] |= 1<<j; a[j] |= 1<<i; }
    return nn;
}

/* record the induced permutation of E(G) if perm is an automorphism */
static void try_perm(void)
{
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            int e1 = ((adj[i]>>j)&1), e2 = ((adj[perm[i]]>>perm[j])&1);
            if (e1 != e2) return;
        }
    if (naut >= MAXAUT) return;
    for (int k = 0; k < nedge; k++)
        aut[naut][k] = pidx[perm[eu[k]]][perm[ev[k]]];
    naut++;
}

static void gen_perm(int d)
{
    if (d == n) { try_perm(); return; }
    for (int v = 0; v < n; v++)
        if (!used[v]) { used[v] = 1; perm[d] = v; gen_perm(d+1); used[v] = 0; }
}

int main(int argc, char **argv)
{
    int MU = (argc > 1) ? atoi(argv[1]) : 2;
    if (MU < 1 || MU > 9) { fprintf(stderr, "usage: %s MU   (1 <= MU <= 9)\n", argv[0]); return 1; }

    char line[256];
    long long total = 0, kept = 0;
    int mult[MAXE], cand[MAXE];

    while (fgets(line, sizeof line, stdin)) {
        size_t L = strlen(line);
        while (L && (line[L-1]=='\n' || line[L-1]=='\r')) line[--L] = 0;
        if (!L || line[0] == '>') continue;
        n = g6_decode(line, adj);
        if (n < 0) { fprintf(stderr, "bad graph6: %s\n", line); return 1; }

        nedge = 0;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) pidx[i][j] = -1;
        for (int i = 0; i < n; i++)
            for (int j = i+1; j < n; j++)
                if ((adj[i]>>j)&1) { eu[nedge]=i; ev[nedge]=j; pidx[i][j]=pidx[j][i]=nedge; nedge++; }
        if (nedge == 0) { printf("%d 0\n", n); kept++; total++; continue; }

        naut = 0; memset(used, 0, sizeof used); gen_perm(0);

        for (int k = 0; k < nedge; k++) mult[k] = 1;
        for (;;) {
            total++;
            /* canonical iff lexicographically least in its Aut(G) orbit */
            int is_min = 1;
            for (int a = 0; a < naut && is_min; a++) {
                for (int k = 0; k < nedge; k++) cand[k] = mult[aut[a][k]];
                for (int k = 0; k < nedge; k++) {
                    if (cand[k] < mult[k]) { is_min = 0; break; }
                    if (cand[k] > mult[k]) break;
                }
            }
            if (is_min) {
                kept++;
                printf("%d %d", n, nedge);
                for (int k = 0; k < nedge; k++)
                    printf(" %d %d %d", eu[k], ev[k], mult[k]);
                putchar('\n');
            }
            int k = nedge - 1;
            while (k >= 0 && mult[k] == MU) mult[k--] = 1;
            if (k < 0) break;
            mult[k]++;
        }
    }
    fprintf(stderr, "gen_multi: %lld assignments, %lld pairwise non-isomorphic\n",
            total, kept);
    return 0;
}
