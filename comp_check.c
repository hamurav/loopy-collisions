/*  loopy_collisions.c   --   PROFILE-ASSUMING VARIANT
 *
 *  ###################################################################
 *  #  THIS PROGRAM ASSUMES AN UNPUBLISHED LEMMA:                     #
 *  #                                                                 #
 *  #      L_G  determines the induced edge-count profile             #
 *  #          #{ B subset V : |B| = k,  e_G(B) = j }                 #
 *  #                                                                 #
 *  #  (iterated vertex deletion).  If that lemma is false, this      #
 *  #  program SILENTLY MISSES COLLISIONS -- it does not error.       #
 *  #  Use loopy_collisions (no -P) for anything you intend to        #
 *  #  publish until the lemma is written down and checked.           #
 *  ###################################################################
 *
 *  Note that the lemma cannot be tested at n <= 10: there, equal L implies
 *  equal Lhat (verified) implies equal polychromate implies equal profile,
 *  so agreement is automatic and carries no information.  Reproducing the
 *  known n <= 9 answers below therefore checks the CODE, not the LEMMA.
 *
 *  Differences from loopy_collisions.c:
 *    * the profile is part of stage 0 and is on by default (--no-profile off);
 *    * the matching-polynomial stage is gone -- after the profile it prunes
 *      nothing at all (28 -> 28 at n=8, 474 -> 474 at n=9);
 *    * --shard i/N shards on the profile key, which is far finer than the
 *      degree sequence, so the shards come out much more evenly sized.
 *
 *  Sieve:  stage 0  m, degree sequence, PROFILE             (~2^n)
 *          stage 1  + #spanning trees, #triangles           (~n^3)
 *          stage 2  + chromatic symmetric function X_G      (2 * 3^n)
 *          stage 3  full L / Lhat / XB fingerprints         (~12 * 3^n)
 *
 *  Measured survival (complete lists):  n=8  0.252% then 0.144%
 *                                       n=9  0.182% then 0.058%
 *
 *  BUILD    cc -O3 -pthread -o loopy_collisions loopy_collisions.c
 *  RUN      geng -q -c 11 27:27 | ./loopy_collisions -q
 *
 *  Validation checkpoints (connected simple graphs):
 *      n <= 7 : no collisions
 *      n = 8  : 8 classes, all pairs, 13..15 edges, L == Lhat == XB
 *      n = 9  : 65 classes, all pairs, 13..23 edges, L == Lhat == XB
 */

#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
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
 * Do NOT add filters that are only known for XB (independence polynomial, the
 * induced-edge-count profile of the polychromate, ...): they could hide exactly
 * the counterexample the search is looking for.
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


/* ============================================================================
 * comp_check.c -- is the set of collisions closed under complementation?
 *
 * If L_G determines L_{Gbar} (Problem 7.1 / Outlook item 2), then
 *      L_G = L_H   ==>   L_{Gbar} = L_{Hbar},
 * so every collision class must map to a collision class under complementation.
 * This program tests that directly on the classes an actual run found.
 *
 * BUILD  cc -O2 -o comp_check comp_check.c -lpthread
 * RUN    ./comp_check < res11/all.txt          (reads the CLASS lines)
 *        ./comp_check < collision_pairs.txt    (also reads "n m g6 g6" lines)
 *        ./comp_check seed < ...               (repeat with a different seed)
 *
 * Complements of connected graphs need not be connected; the "conn" column
 * flags that, because the original search covered connected graphs only and a
 * disconnected complementary class would not have appeared in its output.
 * ==========================================================================*/
static int is_connected(const int *adj,int n){
    if(n==0) return 1;
    int seen=1, frontier=1;
    while(frontier){ int nx=0;
        for(int v=0;v<n;v++) if((frontier>>v)&1) nx|=adj[v];
        nx&=~seen; seen|=nx; frontier=nx; }
    return seen==((1<<n)-1);
}
int main(int argc,char**argv)
{
    unsigned long long seed = (argc>1)?strtoull(argv[1],0,10):12345;
    K = 8;
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


    Work *w = (Work*)malloc(sizeof(Work));
    char line[8192];
    long ncls=0, ok=0, bad=0, bad_conn=0, bad_disc=0;
    printf("%-5s %-4s %-4s %-6s %-6s  %s\n","size","m","mbar","class?","conn?","members (complemented)");
    while (fgets(line,sizeof line,stdin)) {
        char *g6[64]; int ng=0, n=-1, m=-1;
        if (!strncmp(line,"CLASS",5)) {
            char *c = strchr(line,':'); if(!c) continue;
            char *mm = strstr(line,"m="); if(mm) m=atoi(mm+2);
            for(char*t=strtok(c+1," \t\n"); t && ng<64; t=strtok(0," \t\n")) g6[ng++]=t;
        } else if (line[0]=='#'||line[0]=='\n') { continue;
        } else {                                  /* "n m g6 g6 ..." */
            char *t=strtok(line," \t\n"); if(!t||!isdigit((unsigned char)*t)) continue;
            n=atoi(t); t=strtok(0," \t\n"); if(!t) continue; m=atoi(t);
            for(t=strtok(0," \t\n"); t && ng<64; t=strtok(0," \t\n")) g6[ng++]=t;
        }
        if (ng<2) continue;
        int A[64][MAXN]; u64 fL[64],fH,fX; int mm2=0, allconn=1, same=1;
        int nn=0;
        for(int k=0;k<ng;k++){
            int adj[MAXN]; nn=g6_decode(g6[k],adj); if(nn<0){ng=0;break;}
            int full=(1<<nn)-1;
            for(int v=0;v<nn;v++) A[k][v]=(~adj[v])&full&~(1<<v);
            if(!is_connected(A[k],nn)) allconn=0;
            fingerprints(w,A[k],nn,&fL[k],&fH,&fX,&mm2);
        }
        if(ng<2) continue;
        for(int k=1;k<ng;k++) if(fL[k]!=fL[0]) same=0;
        ncls++;
        if(same) ok++; else { bad++; if(allconn) bad_conn++; else bad_disc++; }
        if(!same || ncls<=6){
            printf("%-5d %-4d %-4d %-6s %-6s ",ng,m,mm2,same?"yes":"NO",allconn?"conn":"disc");
            for(int k=0;k<ng;k++) printf(" %s",g6[k]);
            printf("%s\n", same?"":"   <== complements do NOT collide");
        }
    }
    printf("\nclasses tested: %ld\n",ncls);
    printf("  complements still collide : %ld\n",ok);
    printf("  complements do NOT collide: %ld   (of these, %ld have all complements connected)\n",bad,bad_conn);
    if(bad_conn) printf("\n==> L_G does NOT determine L_Gbar: a connected counterexample above.\n");
    else if(bad)  printf("\n==> counterexamples exist but involve disconnected complements;\n"
                         "    still decisive, since L is multiplicative over components.\n");
    else          printf("\n==> no counterexample in this input: collisions are closed under complementation here.\n");
    return 0;
}
