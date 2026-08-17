/*  mg_collisions.c  --  collision search for loopless MULTIGRAPHS
 *
 *  Companion to loopy_collisions.c (which handles simple graphs).  Nothing in
 *  that program is changed or reused at run time; this is a separate tool.
 *
 *  PURPOSE
 *    Conjecture 4.2 is stated for simple graphs.  Example 4.10 shows that it
 *    fails if loops are allowed.  The status for loopless multigraphs -- i.e.
 *    whether parallel edges are an obstruction -- is open.  This program looks
 *    for a counterexample: a pair of non-isomorphic loopless multigraphs with
 *    equal L but different Lhat, or equal XB but different L, etc.
 *
 *  WHY THE ALGORITHM IS UNCHANGED
 *    Every invariant here is computed from the connected-block polynomials
 *    C_G(B;q), and those depend on G only through the edge counts e_G(S) of
 *    induced subgraphs.  Counting edges WITH MULTIPLICITY is the only change
 *    needed; the O(3^n) recursion and the set-partition transforms are
 *    identical.
 *
 *  INPUT   the text output of nauty's multig, i.e. whitespace-separated
 *          n  m  v1 w1 mult1  v2 w2 mult2  ...  (repeated, one graph per record)
 *
 *      geng -q -c 7 | multig -m3 -T | ./mg_collisions
 *
 *          multig assigns multiplicities 1..mu to the edges of each simple
 *          graph produced by geng, and emits pairwise non-isomorphic results.
 *          Hence ANY class of size >= 2 found here is a genuine collision --
 *          no isomorphism testing is needed in this program.
 *
 *  NO SIEVE
 *    Unlike loopy_collisions.c this program computes the full fingerprints for
 *    every input graph.  The cheap filters of the simple-graph sieve would all
 *    have to be re-justified for multigraphs, and at the vertex counts that are
 *    reachable here the sieve buys nothing: 3^7 = 2187.
 *
 *  BUILD   cc -O3 -pthread -o mg_collisions mg_collisions.c
 *  OPTIONS -k K   substitution points (default 8)
 *          -s S   random seed; rerun with a second seed as a check
 *          -q     quiet
 *
 *  MEMORY  about 60 bytes per input graph; ~10^7 graphs is comfortable.
 */

/*  loopy_collisions.c
 *
 *  Exhaustive search for pairs of non-isomorphic connected simple graphs with
 *  equal loopy polynomial L_G, refined loopy polynomial Lhat_G, or Tutte
 *  symmetric function XB_G.  Companion code for
 *
 *      A. Kirillov, G. Nenashev, B. Shapiro, A. Vaintrob,
 *      "The Loopy Polynomial: from Tutte's Universal V-Function to
 *       Bizonotopal Geometry", Section 4.3.
 *
 *  Verified through n = 11 vertices (1 006 700 565 connected graphs).
 *
 *  METHOD
 *    Connected-block polynomials C_G(B;q) for all 2^n vertex subsets by the
 *    standard subset recursion (O(3^n) ring operations), then set-partition
 *    dynamic programming for each invariant.  The expansion of C_G(B;q) in the
 *    basis q^{s-1}(1+q)^r is avoided by the identity
 *
 *        sum_r a_r theta^r  =  C_G(B; theta-1) / (theta-1)^{s-1},
 *
 *    so a substitution u_{s,r} -> sum_k beta_{k,s} theta_k^r needs only K
 *    numerical evaluations.  All arithmetic is in F_p, p = 2^61 - 1, at random
 *    substitution points.  Fingerprint equality is a NECESSARY condition for
 *    equality of the invariants, so no collision is missed; the survivors are
 *    confirmed exactly over Z by exact_check.py.
 *
 *  CHEAP SIEVE  (each stage applied only to the survivors of the previous one;
 *  every filter is determined by L_G AND by XB_G, so nothing that could be a
 *  collision is ever discarded)
 *
 *      stage 0   m, degree sequence, induced edge-count profile      ~2^n
 *      stage 1   + number of spanning trees, number of triangles     ~n^3
 *      stage 2   + chromatic symmetric function X_G = XB_G(-1)      2*3^n
 *      stage 3   full L / Lhat / XB fingerprints                   12*3^n
 *
 *  Measured survival on complete lists:  n=8  0.252% -> 0.144%
 *                                        n=9  0.182% -> 0.058%
 *
 *  BUILD   cc -O3 -pthread -o loopy_collisions loopy_collisions.c
 *          (macOS: use cc/clang; do NOT pass -fopenmp, this code uses pthreads)
 *  RUN     geng -q -c 11 27:27 | ./loopy_collisions -q
 *
 *  OPTIONS
 *      -k K            number of substitution points (default 8)
 *      -s SEED         random seed; rerun with a second seed as a check
 *      --no-profile    drop the profile from stage 0 (slower; for cross-checks)
 *      --shard i/N     process only shard i of N, split on the stage-0 key
 *      -o FILE         dump the fingerprint records
 *      -q              quiet
 *
 *  VALIDATION CHECKPOINTS (connected simple graphs)
 *      n <= 7 : no collisions
 *      n = 8  : 8 classes, all of size 2, 13..15 edges, L == Lhat == XB
 *      n = 9  : 65 classes, all of size 2, 13..23 edges, L == Lhat == XB
 *      n = 10 : 1285 classes (1281 of size 2, 4 of size 4), 14..32 edges
 *      n = 11 : 22499 classes (22491 of size 2, 4 of size 3, 4 of size 4),
 *               13..42 edges
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <unistd.h>

typedef unsigned long long u64;
typedef __int128 u128;

/* 2^61 - 1 (Mersenne).  gcc turns the __int128 reduction into a few ops. */
static const u64 P = 2305843009213693951ULL;
static inline u64 mulm(u64 a, u64 b){ return (u64)(((u128)a*b) % P); }
static inline u64 addm(u64 a, u64 b){ u64 s=a+b; return s>=P ? s-P : s; }
static inline u64 subm(u64 a, u64 b){ return a>=b ? a-b : a+P-b; }
static u64 powm(u64 a, u64 e){ u64 r=1; while(e){ if(e&1) r=mulm(r,a); a=mulm(a,a); e>>=1; } return r; }
static inline u64 invm(u64 a){ return powm(a, P-2); }

#define MAXN   12   /* supports n <= 12; keeps scratch arrays cache-friendly */
#define MAXK   24
#define MAXDEG (MAXN*(MAXN-1)/2 + 2)

static int  K = 8;                       /* number of geometric components   */
static int  shard_i = 0, shard_N = 1;    /* --shard i/N : keep 1/N of the input */
static int  use_profile = 1;             /* --no-profile disables (then use the other program) */
static u64  q0, pi_[MAXN+2];             /* XB   random point                */
static u64  theta[MAXK];
static u64  betaL[MAXK][MAXN+2];         /* diagonal (loopy) substitution    */
static u64  betaH[MAXK][MAXN+2];
static u64  yU, iyU[MAXN+2], zU[MAXN+2], pwU[MAXDEG];         /* generic  (refined) substitution  */
static u64  pw_q0[MAXDEG];               /* (1+q0)^i                         */
static u64  pw_th[MAXK][MAXDEG];         /* (1+(theta_k-1))^i = theta_k^i    */
static u64  iqp[MAXK][MAXN+2];           /* (theta_k-1)^{-(s-1)}             */

/* ------------------------------------------------------------------ graph6 */
static int g6_decode(const char *s, int *adj)
{
    int n = (int)s[0] - 63;
    if (n < 1 || n > MAXN) return -1;
    for (int i = 0; i < n; i++) adj[i] = 0;
    int nbits = n*(n-1)/2, p = 0;
    const char *c = s + 1;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++, p++) {
            int byte = p / 6, bit = 5 - (p % 6);
            int v = (int)c[byte] - 63;
            if (v < 0 || v > 63) return -1;
            if ((v >> bit) & 1) { adj[i] |= 1<<j; adj[j] |= 1<<i; }
        }
    (void)nbits;
    return n;
}

/* shard key: m together with the sorted degree sequence.  Both L_G and XB_G
 * determine these (Cor 2.11(3), Thm 3.10), so two graphs that could possibly
 * collide always receive the same key and hence land in the same shard.       */
static u64 sieve_edge_profile(const int *adj, int n);
static u64 shard_key(const int *adj, int n)
{
    int deg[16], m = 0;
    for (int v = 0; v < n; v++) { deg[v] = __builtin_popcount(adj[v]); m += deg[v]; }
    m /= 2;
    for (int a = 0; a < n; a++)               /* insertion sort */
        for (int b = a+1; b < n; b++)
            if (deg[b] < deg[a]) { int t = deg[a]; deg[a] = deg[b]; deg[b] = t; }
    u64 h = 1469598103934665603ULL;
    h = h*1000003ULL + (u64)m;
    for (int v = 0; v < n; v++) h = h*1000003ULL + (u64)deg[v];
    h = h*1000003ULL + sieve_edge_profile(adj, n);   /* far finer than deg seq */
    h ^= h >> 29; h *= 0xbf58476d1ce4e5b9ULL; h ^= h >> 32;
    return h;
}

/* ------------------------------------------------------------- per-graph DP */
typedef struct { int ec[1<<MAXN]; u64 N[1<<MAXN], C[1<<MAXN], g[1<<MAXN], Z[1<<MAXN]; } Work;

/* -------------------------------------------------- portable parallel-for --
 * pthreads only: no OpenMP, so this builds with Apple clang as well as gcc.
 * Each worker owns one Work scratch buffer and grabs chunks from a counter.
 * Thread count: $NTHREADS, else the number of online CPUs, else 4.          */
typedef void (*par_body)(size_t i, Work *w, void *ctx);
typedef struct {
    size_t n, chunk, next;
    pthread_mutex_t mu;
    par_body body;
    void *ctx;
} ParJob;

static void *par_worker(void *arg)
{
    ParJob *J = (ParJob*)arg;
    Work *w = (Work*)malloc(sizeof(Work));
    if (!w) { fprintf(stderr,"out of memory (Work)\n"); exit(1); }
    for (;;) {
        pthread_mutex_lock(&J->mu);
        size_t lo = J->next; J->next += J->chunk;
        pthread_mutex_unlock(&J->mu);
        if (lo >= J->n) break;
        size_t hi = lo + J->chunk; if (hi > J->n) hi = J->n;
        for (size_t i = lo; i < hi; i++) J->body(i, w, J->ctx);
    }
    free(w);
    return NULL;
}

static int par_nthreads(void)
{
    const char *e = getenv("NTHREADS");
    if (e) { int t = atoi(e); if (t > 0) return t; }
#ifdef _SC_NPROCESSORS_ONLN
    long k = sysconf(_SC_NPROCESSORS_ONLN);
    if (k > 0) return (int)k;
#endif
    return 4;
}

static void parallel_for(size_t n, par_body body, void *ctx, size_t chunk)
{
    if (n == 0) return;
    int T = par_nthreads();
    if ((size_t)T > (n + chunk - 1)/chunk) T = (int)((n + chunk - 1)/chunk);
    if (T < 1) T = 1;
    ParJob J; J.n = n; J.chunk = chunk; J.next = 0; J.body = body; J.ctx = ctx;
    pthread_mutex_init(&J.mu, NULL);
    pthread_t *th = (pthread_t*)malloc(sizeof(pthread_t)*(size_t)T);
    int started = 0;
    for (int i = 1; i < T; i++)
        if (pthread_create(&th[started], NULL, par_worker, &J) == 0) started++;
    par_worker(&J);                       /* the calling thread works too */
    for (int i = 0; i < started; i++) pthread_join(th[i], NULL);
    free(th);
    pthread_mutex_destroy(&J.mu);
}

static void calc_ecount(const int *adj, int n, int *ec)
{
    for (int S = 0; S < (1<<n); S++) ec[S] = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            if ((adj[i]>>j)&1)
                for (int S = 0; S < (1<<n); S++)
                    if (((S>>i)&1) && ((S>>j)&1)) ec[S]++;
}

/* C[S] = connected-spanning-subgraph enumerator of G[S] at the q with (1+q)^i = pw[i] */
static void connenum(Work *w, int n, const u64 *pw)
{
    for (int S = 0; S < (1<<n); S++) w->N[S] = pw[w->ec[S]];
    w->C[0] = 0;
    for (int S = 1; S < (1<<n); S++) {
        int low = S & (-S), rest = S ^ low;
        u64 acc = w->N[S];
        for (int T = (rest-1)&rest; ; T = (T-1)&rest) {
            int B = T | low;
            if (B != S) acc = subm(acc, mulm(w->C[B], w->N[S^B]));
            if (T == 0) break;
        }
        w->C[S] = acc;
    }
}

/* Z[V] = sum over set partitions pi of prod_{B in pi} g[B] */
static u64 combine(Work *w, int n)
{
    w->Z[0] = 1;
    for (int S = 1; S < (1<<n); S++) {
        int low = S & (-S), rest = S ^ low;
        u64 acc = 0;
        for (int T = rest; ; T = (T-1)&rest) {
            int B = T | low;
            acc = addm(acc, mulm(w->g[B], w->Z[S^B]));
            if (T == 0) break;
        }
        w->Z[S] = acc;
    }
    return w->Z[(1<<n)-1];
}

/* fingerprints of L, Lhat and XB, in one sweep over the K theta values */
static void fingerprints(Work *w, const int *adj, int n, u64 *fL, u64 *fH, u64 *fX, int *medges)
{
    calc_ecount(adj, n, w->ec);
    *medges = w->ec[(1<<n)-1];

    static __thread u64 gL[1<<MAXN], gH[1<<MAXN];
    for (int S = 0; S < (1<<n); S++) { gL[S] = 0; gH[S] = 0; }

    for (int k = 0; k < K; k++) {
        connenum(w, n, pw_th[k]);
        for (int S = 1; S < (1<<n); S++) {
            int s = __builtin_popcount(S);
            u64 base = mulm(w->C[S], iqp[k][s]);     /* C(B;theta-1)/(theta-1)^{s-1} */
            gL[S] = addm(gL[S], mulm(betaL[k][s], base));
            gH[S] = addm(gH[S], mulm(betaH[k][s], base));
        }
    }
    memcpy(w->g, gL, sizeof(u64)*(1<<n)); *fL = combine(w, n);
    memcpy(w->g, gH, sizeof(u64)*(1<<n)); *fH = combine(w, n);

    connenum(w, n, pw_q0);
    w->g[0] = 0;
    for (int S = 1; S < (1<<n); S++)
        w->g[S] = mulm(w->C[S], pi_[__builtin_popcount(S)]);
    *fX = combine(w, n);
}

/* ============================ CHEAP SIEVE ==================================
 * Every filter below is implied by BOTH L_G and XB_G, so a graph pruned here
 * cannot belong to a collision class of either invariant:
 *      m                     Cor 2.11(3);  top q-degree of XB
 *      degree sequence       Thm 3.10;     classical on the U / polychromate side
 *      #spanning trees       = T_G(1,1);   both determine T_G (Cor 2.4)
 *      #triangles            Whitney, from the chromatic polynomial = T_G spec.
 *      matching polynomial   Cor 3.13;     [NW] on the U side
 *      X_G                   Thm 3.11;     XB_G(-1)
 * The induced edge-count profile is admissible because L_G determines it
 * (iterated vertex deletion).  Do NOT add filters that are known only for XB
 * (e.g. the independence polynomial): those could hide exactly the
 * counterexample the search is looking for.
 * ========================================================================== */
static int cmp_int_(const void *a, const void *b){ return *(const int*)a - *(const int*)b; }

static u64 sieve_trees(const int *adj, int n)          /* Matrix-Tree, mod P */
{
    if (n == 1) return 1;
    static __thread u64 L[MAXN][MAXN];
    for (int i = 0; i < n-1; i++)
        for (int j = 0; j < n-1; j++)
            L[i][j] = (i==j) ? (u64)__builtin_popcount(adj[i])
                             : (((adj[i]>>j)&1) ? P-1 : 0);
    u64 det = 1;
    for (int c = 0; c < n-1; c++) {
        int piv = -1;
        for (int r = c; r < n-1; r++) if (L[r][c]) { piv = r; break; }
        if (piv < 0) return 0;
        if (piv != c) { for (int j=c;j<n-1;j++){u64 t=L[c][j];L[c][j]=L[piv][j];L[piv][j]=t;} det=subm(0,det); }
        det = mulm(det, L[c][c]);
        u64 ic = invm(L[c][c]);
        for (int r = c+1; r < n-1; r++) if (L[r][c]) {
            u64 f = mulm(L[r][c], ic);
            for (int j = c; j < n-1; j++) L[r][j] = subm(L[r][j], mulm(f, L[c][j]));
        }
    }
    return det;
}

static u64 sieve_matching(const int *adj, int n)       /* matching polynomial */
{
    static __thread u64 f[1<<MAXN][MAXN/2+2];
    int Kc = n/2 + 1;
    for (int k = 0; k < Kc; k++) f[0][k] = 0;
    f[0][0] = 1;
    for (int S = 1; S < (1<<n); S++) {
        int v = __builtin_ctz(S), rest = S ^ (1<<v);
        for (int k = 0; k < Kc; k++) f[S][k] = f[rest][k];
        int nb = adj[v] & rest;
        while (nb) {
            int u = __builtin_ctz(nb); nb &= nb-1;
            int T = rest ^ (1<<u);
            for (int k = 1; k < Kc; k++) f[S][k] = addm(f[S][k], f[T][k-1]);
        }
    }
    u64 h = 1469598103934665603ULL;
    for (int k = 0; k < Kc; k++) h = h*1000003ULL + f[(1<<n)-1][k];
    return h;
}

static u64 sieve_XG(Work *w, int n)                    /* X_G  =  XB_G(-1)   */
{
    for (int S = 0; S < (1<<n); S++) w->N[S] = w->ec[S] ? 0 : 1;   /* (1+q)^e at q=-1 */
    w->C[0] = 0;
    for (int S = 1; S < (1<<n); S++) {
        int low = S&(-S), rest = S^low; u64 acc = w->N[S];
        for (int T = (rest-1)&rest; ; T = (T-1)&rest) {
            int B = T|low;
            if (B != S) acc = subm(acc, mulm(w->C[B], w->N[S^B]));
            if (T == 0) break;
        }
        w->C[S] = acc;
    }
    w->g[0] = 0;
    for (int S = 1; S < (1<<n); S++)
        w->g[S] = mulm(w->C[S], pi_[__builtin_popcount(S)]);
    return combine(w, n);
}

/* induced edge-count profile  #{B : |B|=k, e_G(B)=j}.
 * Legitimate ONLY given the theorem that L_G determines it (iterated vertex
 * deletion); it is classical on the XB / polychromate side.  Enabled by -P.
 * Cost: 2^n popcounts + 2^n histogram increments -- essentially free.        */
static u64 prof_rnd[MAXN+1][MAXN*(MAXN-1)/2+2];   /* filled in main() */
static u64 sieve_edge_profile(const int *adj, int n)
{
    /* Order-independent multiset hash of { (|S|, e_G(S)) : S subset V }, which
     * is exactly the profile.  No histogram array, no memset: one table lookup
     * and one add per subset.  A hash clash can only merge two profiles, i.e.
     * create extra work downstream -- it can never split a genuine class.     */
    static __thread int ec[1<<MAXN];
    ec[0] = 0;
    u64 h = prof_rnd[0][0];
    for (int S = 1; S < (1<<n); S++) {
        int v = __builtin_ctz(S), rest = S ^ (1<<v);
        ec[S] = ec[rest] + __builtin_popcount(adj[v] & rest);
        h += prof_rnd[__builtin_popcount(S)][ec[S]];
    }
    return h;
}

/* stage 0: m, degree sequence, PROFILE.
 * Ordered by value per microsecond: on the complete n=9 list the profile alone
 * takes 99.8% -> 9.2% for ~5.6 us, whereas #spanning trees + #triangles take
 * 99.8% -> 51% for ~5.3 us.  So the profile goes first and the (still useful)
 * spanning-tree count is only paid for by the ~9% that survive it.            */
static u64 sieve_stage0(const int *adj, int n)
{
    int deg[MAXN], m = 0;
    for (int v = 0; v < n; v++) { deg[v] = __builtin_popcount(adj[v]); m += deg[v]; }
    m /= 2;
    qsort(deg, n, sizeof(int), cmp_int_);
    u64 h = 1469598103934665603ULL;
    h = h*1000003ULL + (u64)m;
    for (int v = 0; v < n; v++) h = h*1000003ULL + (u64)deg[v];
    if (use_profile) h = h*1000003ULL + sieve_edge_profile(adj, n);
    return h;
}

/* stage 1: + #spanning trees (Matrix-Tree) and #triangles  (~n^3 operations) */
static u64 sieve_trees_tri(const int *adj, int n, u64 prev)
{
    int tri = 0;
    for (int a = 0; a < n; a++)
        for (int b = a+1; b < n; b++)
            if ((adj[a]>>b)&1)
                tri += __builtin_popcount(adj[a] & adj[b] & ~((1<<(b+1))-1));
    return (prev*1000003ULL + sieve_trees(adj, n))*1000003ULL + (u64)tri;
}
/* stage 1: matching polynomial            (~2^n * n * n/2 operations)        */
/* stage 2: chromatic symmetric function   (~2 * 3^n operations)              */
static u64 sieve_stage(Work *w, const int *adj, int n, int stage, u64 prev)
{
    if (stage == 0) return sieve_stage0(adj, n);            /* m, deg, profile */
    if (stage == 1) return sieve_trees_tri(adj, n, prev);   /* trees, triangles */
    calc_ecount(adj, n, w->ec);                             /* X_G              */
    return prev*1000003ULL + sieve_XG(w, n);
}

static int cmpu64(const void*a,const void*b){u64 x=*(const u64*)a,y=*(const u64*)b;return x<y?-1:(x>y);}
typedef struct { u64 k; unsigned i; } KI;
static int cmpKI(const void *a, const void *b){
    u64 x = ((const KI*)a)->k, y = ((const KI*)b)->k;
    return x < y ? -1 : (x > y);
}

/* contexts for the two parallel loops */
typedef struct { char *g6s; int G6W, n, stage; KI *ki; } SieveCtx;
static void sieve_body(size_t t, Work *w, void *vctx)
{
    SieveCtx *c = (SieveCtx*)vctx;
    int adj[MAXN];
    g6_decode(c->g6s + (size_t)c->ki[t].i * c->G6W, adj);
    c->ki[t].k = sieve_stage(w, adj, c->n, c->stage, c->ki[t].k);
}

/* ------------------------------------------------------------------ records */
typedef struct { u64 fL, fH, fX; unsigned idx; int m; } Rec;
typedef struct { char *g6s; int G6W, n; unsigned *keep; Rec *recs; } FpCtx;
static Rec  *recs;
static char *g6s;          /* fixed-width graph6 storage */
static int   G6W = 24;

static int cmp_by(const void *a, const void *b, int which)
{
    const Rec *x = a, *y = b;
    u64 kx = which==0 ? x->fL : which==1 ? x->fH : x->fX;
    u64 ky = which==0 ? y->fL : which==1 ? y->fH : y->fX;
    return kx < ky ? -1 : kx > ky;
}
static int cmpL(const void*a,const void*b){return cmp_by(a,b,0);}
static int cmpH(const void*a,const void*b){return cmp_by(a,b,1);}
static int cmpX(const void*a,const void*b){return cmp_by(a,b,2);}
static int cmpU64(const void*a,const void*b){u64 x=*(const u64*)a,y=*(const u64*)b;return x<y?-1:(x>y);}

static void fp_body(size_t t, Work *w, void *vctx)
{
    FpCtx *c = (FpCtx*)vctx;
    int adj[MAXN];
    size_t i = c->keep ? (size_t)c->keep[t] : t;
    int nn = g6_decode(c->g6s + i*c->G6W, adj);
    if (nn != c->n) { fprintf(stderr,"bad graph6 at line %zu\n", i); exit(1); }
    u64 fL,fH,fX; int m;
    fingerprints(w, adj, c->n, &fL, &fH, &fX, &m);
    c->recs[t].fL=fL; c->recs[t].fH=fH; c->recs[t].fX=fX;
    c->recs[t].idx=(unsigned)i; c->recs[t].m=m;
}


/* ------------------------------------------------------------ multigraph IO */
#define MAXPAIR (MAXN*(MAXN-1)/2)
static int NV;                       /* vertices, from the first record       */
static unsigned char *mults;         /* NV*(NV-1)/2 multiplicities per graph  */
static size_t ngraphs, capg;
static int STRIDE = MAXPAIR;   /* set to NV*(NV-1)/2 once NV is known */

static int pair_index(int i,int j)   /* i<j */
{
    int k=0;
    for(int a=0;a<i;a++) k += NV-1-a;
    return k + (j-i-1);
}

/* read one whitespace-separated integer; returns 0 at EOF */
static int rdint(long *out)
{
    int c, neg=0; long v=0;
    do { c=getchar(); if(c==EOF) return 0; } while(c==' '||c=='\n'||c=='\t'||c=='\r');
    if(c=='-'){neg=1;c=getchar();}
    if(c<'0'||c>'9') return 0;
    while(c>='0'&&c<='9'){ v=v*10+(c-'0'); c=getchar(); }
    *out = neg? -v : v;
    return 1;
}

/* e_G(S) with multiplicity, for every subset S */
static void calc_ecount_multi(const unsigned char *mu, int n, int *ec)
{
    for(int S=0;S<(1<<n);S++) ec[S]=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
        int k = mu[pair_index(i,j)];
        if(!k) continue;
        for(int S=0;S<(1<<n);S++)
            if(((S>>i)&1)&&((S>>j)&1)) ec[S]+=k;
    }
}

/* fingerprints, identical to the simple-graph code except for the edge counts */
static void fingerprints_multi(Work *w, const unsigned char *mu, int n,
                               u64 *fL, u64 *fH, u64 *fX, u64 *fU, int *medges)
{
    calc_ecount_multi(mu, n, w->ec);
    *medges = w->ec[(1<<n)-1];

    static __thread u64 gL[1<<MAXN], gH[1<<MAXN];
    for (int S = 0; S < (1<<n); S++) { gL[S]=0; gH[S]=0; }
    for (int k = 0; k < K; k++) {
        connenum(w, n, pw_th[k]);
        for (int S = 1; S < (1<<n); S++) {
            int s = __builtin_popcount(S);
            u64 base = mulm(w->C[S], iqp[k][s]);
            gL[S] = addm(gL[S], mulm(betaL[k][s], base));
            gH[S] = addm(gH[S], mulm(betaH[k][s], base));
        }
    }
    memcpy(w->g, gL, sizeof(u64)*(1<<n)); *fL = combine(w, n);
    memcpy(w->g, gH, sizeof(u64)*(1<<n)); *fH = combine(w, n);

    connenum(w, n, pw_q0);
    w->g[0]=0;
    for (int S=1;S<(1<<n);S++)
        w->g[S] = mulm(w->C[S], pi_[__builtin_popcount(S)]);
    *fX = combine(w, n);

    /* ordinary U-polynomial: u_{s,r} -> z_s y^r is the K=1 case of the same
       evaluation trick, with theta = y and beta_s = z_s.                    */
    connenum(w, n, pwU);
    w->g[0]=0;
    for (int S=1;S<(1<<n);S++){
        int s=__builtin_popcount(S);
        w->g[S] = mulm(mulm(w->C[S], iyU[s]), zU[s]);
    }
    *fU = combine(w, n);
}

typedef struct { u64 fL,fH,fX,fU; unsigned idx; int m; } MRec;
static MRec *mrecs;
static int cmpML(const void*a,const void*b){
    const MRec*x=a,*y=b; return x->fL<y->fL?-1:(x->fL>y->fL);
}
static int cmpMH(const void*a,const void*b){
    const MRec*x=a,*y=b; return x->fH<y->fH?-1:(x->fH>y->fH);
}
static int cmpMX(const void*a,const void*b){
    const MRec*x=a,*y=b; return x->fX<y->fX?-1:(x->fX>y->fX);
}
static int cmpMU(const void*a,const void*b){
    const MRec*x=a,*y=b; return x->fU<y->fU?-1:(x->fU>y->fU);
}

static void print_multi(size_t i)
{
    const unsigned char *mu = mults + i*(size_t)STRIDE;
    printf("[");
    int first=1;
    for(int a=0;a<NV;a++) for(int b=a+1;b<NV;b++){
        int k=mu[pair_index(a,b)];
        if(k){ printf("%s%d-%d:%d", first?"":" ", a,b,k); first=0; }
    }
    printf("]");
}

int main(int argc,char**argv)
{
    unsigned long long seed=12345; int verbose=1;
    int only_m = -1;                 /* -m M : keep only graphs with M edges */
    K=8;
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-k")&&i+1<argc) K=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-s")&&i+1<argc) seed=strtoull(argv[++i],0,10);
        else if(!strcmp(argv[i],"-m")&&i+1<argc) only_m=atoi(argv[++i]);
        else if(!strcmp(argv[i],"-q")) verbose=0;
        else { fprintf(stderr,"usage: %s [-k K] [-s SEED] [-m EDGES] [-q] < multig-T-output\n",argv[0]); return 1; }
    }
    /* ---- random parameters ---- */
    u64 st = seed*6364136223846793005ULL + 1442695040888963407ULL;
    #define RND() (st = st*6364136223846793005ULL + 1442695040888963407ULL, (st>>3)%P)
    q0 = RND();
    for (int s = 0; s <= MAXN+1; s++) pi_[s] = RND();
    for (int k = 0; k <= MAXN; k++)
        for (int j = 0; j <= MAXN*(MAXN-1)/2+1; j++) prof_rnd[k][j] = RND();
    pw_q0[0] = 1;
    for (int i = 1; i < MAXDEG; i++) pw_q0[i] = mulm(pw_q0[i-1], addm(1,q0));
    for (int k = 0; k < K; k++) {
        u64 th = RND(); if (th <= 1) th += 11;
        theta[k] = th;
        u64 q = subm(th,1), iq = invm(q);
        /* (1+q)^i = theta^i */
        pw_th[k][0] = 1;
        for (int i = 1; i < MAXDEG; i++) pw_th[k][i] = mulm(pw_th[k][i-1], th);
        u64 alpha = RND();
        for (int s = 0; s <= MAXN+1; s++) {
            iqp[k][s]   = powm(iq, s>0 ? s-1 : 0);
            betaL[k][s] = mulm(alpha, powm(th, s>0 ? s-1 : 0));   /* diagonal   */
            betaH[k][s] = RND();                                  /* generic    */
        }
    }


    yU = RND(); if (yU<=1) yU += 7;
    { u64 iy = invm(subm(yU,1));
      pwU[0]=1; for(int i=1;i<MAXDEG;i++) pwU[i]=mulm(pwU[i-1], yU);
      for(int s=0;s<=MAXN+1;s++){ iyU[s]=powm(iy, s>0?s-1:0); zU[s]=RND(); } }

    /* ---------------------------------------------------------- read input */
    capg = 1u<<16; mults = NULL; ngraphs = 0;
    long n_, m_, v_, w_, k_;
    while (rdint(&n_)) {
        if (!rdint(&m_)) break;
        if (NV==0) { NV=(int)n_;
                     if (NV>MAXN){fprintf(stderr,"n=%d exceeds MAXN=%d\n",NV,MAXN);return 1;}
                     STRIDE = NV*(NV-1)/2;
                     mults = malloc(capg*(size_t)STRIDE); }
        if ((int)n_!=NV) { fprintf(stderr,"mixed vertex counts (%ld vs %d)\n",n_,NV); return 1; }
        if (NV>MAXN)    { fprintf(stderr,"n=%d exceeds MAXN=%d\n",NV,MAXN); return 1; }
        if (ngraphs==capg) { capg*=2; mults=realloc(mults,capg*(size_t)STRIDE);
                              if(!mults){fprintf(stderr,"out of memory at %zu graphs\n",ngraphs);return 1;} }
        unsigned char *mu = mults + ngraphs*(size_t)STRIDE;
        memset(mu,0,(size_t)STRIDE);
        long medge=0;
        for (long t=0;t<m_;t++) {
            if(!rdint(&v_)||!rdint(&w_)||!rdint(&k_)){fprintf(stderr,"short record\n");return 1;}
            if(v_==w_){fprintf(stderr,"loops are not allowed here\n");return 1;}
            int a=(int)(v_<w_?v_:w_), b=(int)(v_<w_?w_:v_);
            mu[pair_index(a,b)] = (unsigned char)k_;
            medge += k_;
        }
        if (only_m >= 0 && medge != only_m) continue;   /* wrong bucket */
        ngraphs++;
    }
    if(!ngraphs){fprintf(stderr,"no graphs read\n");return 1;}
    if(verbose) fprintf(stderr,"read %zu multigraphs on n=%d vertices (K=%d), %.2f GB\n",
                        ngraphs,NV,K,
                        (ngraphs*(double)(STRIDE+sizeof(MRec)))/1e9);

    /* ------------------------------------------------------ fingerprints */
    mrecs = malloc(sizeof(MRec)*ngraphs);
    Work *w = (Work*)malloc(sizeof(Work));
    for (size_t i=0;i<ngraphs;i++){
        int m;
        fingerprints_multi(w, mults+i*(size_t)STRIDE, NV,
                           &mrecs[i].fL,&mrecs[i].fH,&mrecs[i].fX,&mrecs[i].fU,&m);
        mrecs[i].idx=(unsigned)i; mrecs[i].m=m;
    }

    /* ------------------------------------------------------------ group */
    const char *names[4]={"L","Lhat","XB","U"};
    int (*cmps[4])(const void*,const void*)={cmpML,cmpMH,cmpMX,cmpMU};
    u64 *sig[4]={0,0,0,0}; size_t nsig[4]={0,0,0,0};
    for(int which=0;which<4;which++){
        qsort(mrecs,ngraphs,sizeof(MRec),cmps[which]);
        size_t classes=0,pairs=0; int mmin=1<<30,mmax=-1;
        sig[which]=malloc(sizeof(u64)*(ngraphs/2+2));
        for(size_t i=0;i<ngraphs;){
            u64 key = which==0?mrecs[i].fL:which==1?mrecs[i].fH:
                      which==2?mrecs[i].fX:mrecs[i].fU;
            size_t j=i;
            while(j<ngraphs &&
                  (which==0?mrecs[j].fL:which==1?mrecs[j].fH:
                   which==2?mrecs[j].fX:mrecs[j].fU)==key) j++;
            if(j-i>1){
                classes++; pairs += (j-i)*(j-i-1)/2;
                u64 sg=1469598103934665603ULL;
                unsigned mem[64]; size_t nm=0;
                for(size_t t=i;t<j&&nm<64;t++) mem[nm++]=mrecs[t].idx;
                for(size_t a=0;a<nm;a++) for(size_t b=a+1;b<nm;b++)
                    if(mem[b]<mem[a]){unsigned x=mem[a];mem[a]=mem[b];mem[b]=x;}
                for(size_t a=0;a<nm;a++) sg=sg*1000003ULL+mem[a];
                sig[which][nsig[which]++]=sg;
                printf("CLASS %s size=%zu m=%d :",names[which],j-i,mrecs[i].m);
                for(size_t t=i;t<j;t++){ printf(" "); print_multi(mrecs[t].idx);
                    if(mrecs[t].m<mmin)mmin=mrecs[t].m;
                    if(mrecs[t].m>mmax)mmax=mrecs[t].m; }
                printf("\n");
            }
            i=j;
        }
        printf("SUMMARY %s n=%d graphs=%zu classes=%zu pairs=%zu edges=[%d,%d]\n",
               names[which],NV,ngraphs,classes,pairs,
               mmax<0?0:mmin, mmax<0?0:mmax);
    }
    /* ---- does one partition refine another?  (class containment, not just
            class counts)                                                   */
    {
        struct { int a,b; const char *na,*nb; } tests[] = {
            {0,2,"L","XB"}, {2,0,"XB","L"}, {0,1,"L","Lhat"}, {1,0,"Lhat","L"},
            {2,3,"XB","U"}, {3,2,"U","XB"}, {1,2,"Lhat","XB"}
        };
        for (unsigned tt=0; tt<sizeof tests/sizeof tests[0]; tt++) {
            int A=tests[tt].a, B=tests[tt].b;
            qsort(mrecs,ngraphs,sizeof(MRec),cmps[A]);
            int refines=1;
            for(size_t i=0;i<ngraphs && refines;){
                u64 ka = A==0?mrecs[i].fL:A==1?mrecs[i].fH:A==2?mrecs[i].fX:mrecs[i].fU;
                u64 kb = B==0?mrecs[i].fL:B==1?mrecs[i].fH:B==2?mrecs[i].fX:mrecs[i].fU;
                size_t j=i;
                while(j<ngraphs){
                    u64 a2 = A==0?mrecs[j].fL:A==1?mrecs[j].fH:A==2?mrecs[j].fX:mrecs[j].fU;
                    if(a2!=ka) break;
                    u64 b2 = B==0?mrecs[j].fL:B==1?mrecs[j].fH:B==2?mrecs[j].fX:mrecs[j].fU;
                    if(b2!=kb){ refines=0; break; }
                    j++;
                }
                i=j;
            }
            printf("REFINE   %-4s-classes refine %-4s-classes : %s\n",
                   tests[tt].na, tests[tt].nb, refines?"YES":"NO");
        }
    }

    for(int w2=0;w2<4;w2++) qsort(sig[w2],nsig[w2],8,cmpu64);
    #define SAME(a,b) (nsig[a]==nsig[b] && !memcmp(sig[a],sig[b],8*nsig[a]))
    int eqLH=SAME(0,1), eqXU=SAME(2,3), eqHX=SAME(1,2);
    printf("VERDICT  L-classes    == Lhat-classes: %s\n", eqLH?"YES":"NO  <-- COUNTEREXAMPLE");
    printf("VERDICT  XB-classes   == U-classes   : %s\n", eqXU?"YES":"NO  <-- COUNTEREXAMPLE");
    printf("VERDICT  Lhat-classes == XB-classes  : %s\n", eqHX?"YES":"NO  <-- COUNTEREXAMPLE");
    if(eqLH&&eqHX&&eqXU)
        printf("\nNo counterexample among these multigraphs: all four invariants agree.\n");
    return 0;
}
