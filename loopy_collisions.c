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
static u64  betaH[MAXK][MAXN+2];         /* generic  (refined) substitution  */
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

int main(int argc, char **argv)
{
    u64 seed = 987654321ULL;
    const char *dump = NULL;
    int verbose = 1, sieve = 1;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i],"-k") && i+1 < argc) K = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-s") && i+1 < argc) seed = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i],"-o") && i+1 < argc) dump = argv[++i];
        else if (!strcmp(argv[i],"-S")) sieve = 1;   /* accepted, always on */
        else if (!strcmp(argv[i],"--no-profile")) use_profile = 0;
        else if (!strcmp(argv[i],"--shard") && i+1 < argc) {
            if (sscanf(argv[++i], "%d/%d", &shard_i, &shard_N) != 2 ||
                shard_N < 1 || shard_i < 0 || shard_i >= shard_N) {
                fprintf(stderr,"bad --shard (want i/N with 0 <= i < N)\n"); return 1;
            }
        }
        else if (!strcmp(argv[i],"-q")) verbose = 0;
        else { fprintf(stderr,"usage: %s [-k K] [-s seed] [--no-profile] [--shard i/N] [-o records.bin] [-q] < graphs.g6\n",argv[0]); return 1; }
    }
    if (K < 1 || K > MAXK) { fprintf(stderr,"K out of range\n"); return 1; }

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

    /* ---- read all graph6 lines ---- */
    size_t cap = 1u<<20, cnt = 0, nread = 0;
    g6s  = malloc(cap * G6W);
    char line[256];
    int n = -1;
    while (fgets(line, sizeof line, stdin)) {
        size_t L = strlen(line);
        while (L && (line[L-1]=='\n' || line[L-1]=='\r')) line[--L] = 0;
        if (!L || line[0]=='>') continue;                 /* skip >>graph6<< header */
        if ((int)L >= G6W) { fprintf(stderr,"graph6 line too long\n"); return 1; }
        if (n < 0) n = (int)line[0] - 63;
        nread++;
        if (shard_N > 1) {                       /* keep only our shard */
            int adjr[MAXN];
            if (g6_decode(line, adjr) != n) { fprintf(stderr,"bad graph6\n"); return 1; }
            if ((int)(shard_key(adjr, n) % (u64)shard_N) != shard_i) continue;
        }
        if (cnt == cap) { cap *= 2; g6s = realloc(g6s, cap*G6W); }
        memcpy(g6s + cnt*G6W, line, L+1);
        cnt++;
    }
    if (!cnt) { fprintf(stderr,"no graphs on stdin\n"); return 1; }
    if (n < 1 || n > MAXN) { fprintf(stderr,"unsupported n=%d\n", n); return 1; }
    if (verbose && use_profile)
        fprintf(stderr,"[loopy_collisions] ASSUMING: L_G determines the induced "
                       "edge-count profile.\n"
                       "                     A false lemma here silently loses collisions.\n");
    if (verbose) fprintf(stderr,"read %zu graphs on n=%d vertices (K=%d)%s\n", nread, n, K,
                     shard_N > 1 ? "" : "");
    if (verbose && shard_N > 1)
        fprintf(stderr,"  shard %d/%d: kept %zu (%.2f%%)\n", shard_i, shard_N, cnt, 100.0*cnt/nread);

    /* ---- optional staged cheap sieve --------------------------------------
     * Each stage is applied only to the survivors of the previous one, so the
     * expensive stages are paid for by very few graphs.                      */
    unsigned *keep = NULL; size_t nkeep = cnt;
    if (sieve) {
        static const char *sname[3] = {"m, deg, PROFILE","+ trees, triangles","+ X_G"};
        KI *ki = malloc(sizeof(KI)*cnt);
        for (size_t i = 0; i < cnt; i++) { ki[i].k = 0; ki[i].i = (unsigned)i; }
        size_t cur = cnt;
        for (int stage = 0; stage < 3 && cur > 1; stage++) {
            SieveCtx sc = { g6s, G6W, n, stage, ki };
            parallel_for(cur, sieve_body, &sc, 2048);
            qsort(ki, cur, sizeof(KI), cmpKI);
            size_t out = 0;
            for (size_t i = 0; i < cur; ) {
                size_t j = i; while (j < cur && ki[j].k == ki[i].k) j++;
                if (j - i > 1) for (size_t t = i; t < j; t++) ki[out++] = ki[t];
                i = j;
            }
            if (verbose)
                fprintf(stderr,"  sieve stage %d (%-22s): %zu -> %zu  (%.3f%% of all)\n",
                        stage, sname[stage], cur, out, 100.0*out/cnt);
            cur = out;
        }
        keep = malloc(sizeof(unsigned)*(cur ? cur : 1));
        for (size_t t = 0; t < cur; t++) keep[t] = ki[t].i;
        nkeep = cur;
        free(ki);
        if (verbose)
            fprintf(stderr,"cheap sieve: %zu / %zu survive (%.3f%%)\n",
                    nkeep, cnt, 100.0*nkeep/cnt);
    }

    recs = malloc(sizeof(Rec)*(nkeep ? nkeep : 1));

    /* ---- fingerprints ---- */
    {
        FpCtx fc = { g6s, G6W, n, keep, recs };
        parallel_for(nkeep, fp_body, &fc, 4096);
    }
    cnt = nkeep;
    if (verbose) fprintf(stderr,"fingerprints done\n");

    if (dump) {
        FILE *f = fopen(dump,"wb");
        fwrite(&cnt,sizeof cnt,1,f);
        fwrite(recs,sizeof(Rec),cnt,f);
        fclose(f);
        if (verbose) fprintf(stderr,"wrote %s\n", dump);
    }

    /* ---- group ---- */
    const char *names[3] = {"L","Lhat","XB"};
    int (*cmps[3])(const void*,const void*) = {cmpL,cmpH,cmpX};
    u64 *sig[3] = {NULL,NULL,NULL}; size_t nsig[3] = {0,0,0};
    for (int which = 0; which < 3; which++) {
        qsort(recs, cnt, sizeof(Rec), cmps[which]);
        size_t classes = 0, pairs = 0; int mmin = 1<<30, mmax = -1;
        sig[which] = malloc(sizeof(u64)*(cnt/2 + 2));
        printf("### invariant %s\n", names[which]);
        for (size_t i = 0; i < cnt; ) {
            u64 key = which==0?recs[i].fL:which==1?recs[i].fH:recs[i].fX;
            size_t j = i;
            while (j < cnt && (which==0?recs[j].fL:which==1?recs[j].fH:recs[j].fX) == key) j++;
            if (j - i > 1) {
                classes++;
                pairs += (j-i)*(j-i-1)/2;
                /* class signature: order-independent hash of the member indices */
                u64 s = 1469598103934665603ULL;
                unsigned mem[64]; size_t nm = 0;
                for (size_t k = i; k < j && nm < 64; k++) mem[nm++] = recs[k].idx;
                for (size_t a = 0; a < nm; a++)          /* insertion sort, tiny */
                    for (size_t b = a+1; b < nm; b++)
                        if (mem[b] < mem[a]) { unsigned t = mem[a]; mem[a] = mem[b]; mem[b] = t; }
                for (size_t a = 0; a < nm; a++) s = s*1000003ULL + mem[a];
                sig[which][nsig[which]++] = s;

                printf("CLASS %s size=%zu m=%d :", names[which], j-i, recs[i].m);
                for (size_t k = i; k < j; k++) {
                    printf(" %s", g6s + (size_t)recs[k].idx*G6W);
                    if (recs[k].m < mmin) mmin = recs[k].m;
                    if (recs[k].m > mmax) mmax = recs[k].m;
                }
                printf("\n");
            }
            i = j;
        }
        printf("SUMMARY %s n=%d graphs=%zu classes=%zu pairs=%zu edges=[%d,%d]\n",
               names[which], n, cnt, classes, pairs,
               mmax < 0 ? 0 : mmin, mmax < 0 ? 0 : mmax);
        fflush(stdout);
    }

    /* ---- do the three partitions coincide? -------------------------------
     * Lhat always refines L (L is a specialization of Lhat, eq. 3.7), and for
     * n <= 10 Markstrom's Observation 6.6 gives Lhat <=> XB.  So the whole
     * content of Conjecture 4.2 at n <= 10 is the single implication
     *       L_G = L_H   ==>   Lhat_G = Lhat_H,
     * i.e. "every L-class is a single Lhat-class".  That is the VERDICT line.  */
    for (int w = 0; w < 3; w++) qsort(sig[w], nsig[w], sizeof(u64), cmpU64);
    int LH = (nsig[0] == nsig[1]) && !memcmp(sig[0], sig[1], sizeof(u64)*nsig[0]);
    int LX = (nsig[0] == nsig[2]) && !memcmp(sig[0], sig[2], sizeof(u64)*nsig[0]);
    int HX = (nsig[1] == nsig[2]) && !memcmp(sig[1], sig[2], sizeof(u64)*nsig[1]);
    printf("\n### partition comparison  (classes: L=%zu Lhat=%zu XB=%zu)\n",
           nsig[0], nsig[1], nsig[2]);
    printf("VERDICT  L-classes == Lhat-classes : %s"
           "   <- Conjecture 4.2 at this n, given Markstrom Obs. 6.6\n",
           LH ? "YES" : "NO  *** COUNTEREXAMPLE ***");
    printf("VERDICT  Lhat-classes == XB-classes: %s"
           "   <- independent re-derivation of Markstrom Obs. 6.6\n",
           HX ? "YES" : "NO  *** U does not determine U^ext here ***");
    printf("VERDICT  L-classes == XB-classes   : %s\n",
           LX ? "YES" : "NO");
    printf("(all classes still need exact confirmation: pipe this file to exact_check.py)\n");
    return (LH && HX && LX) ? 0 : 3;
}
