/*  gen_graphs.c  --  fallback generator of all connected simple graphs on n
 *  vertices up to isomorphism, written to stdout in graph6.
 *
 *  USE nauty's  geng -q -c n  IF YOU HAVE IT: it is faster, battle-tested, and
 *  lets you split the work by edge count.  This program exists only so that the
 *  pipeline is self-contained.
 *
 *  Method: extend every graph on n-1 vertices by one new vertex with every
 *  possible neighbourhood, then de-duplicate by a strong composite isomorphism
 *  invariant (1-WL colours; closed-walk counts (A^k)_vv, k<=4; pair walk counts
 *  (A^k)_uv, k<=4).  The invariant is *not* proved complete, so the program
 *  checks its output against the known counts A000088 (all graphs) and A001349
 *  (connected).  Because the key is an isomorphism invariant, the number of
 *  buckets can never EXCEED the number of isomorphism classes; therefore
 *  agreement with A000088 certifies that no two classes were merged, i.e. the
 *  de-duplication is exact.  The program aborts if the counts disagree.
 *
 *  BUILD  gcc -O3 -march=native -o gen_graphs gen_graphs.c
 *  RUN    ./gen_graphs 9 > c9.g6
 *         ./gen_graphs 10 > c10.g6      (~25 min, ~1 GB RAM)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
#define MAXN 10

static const u64 A000088[] = {1,1,2,4,11,34,156,1044,12346,274668,12005168};
static const u64 A001349[] = {1,1,1,2, 6,21,112, 853,11117,261080,11716571};

static int PIDX[MAXN][MAXN];
static void init_pidx(int n){int c=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){PIDX[i][j]=c;PIDX[j][i]=c;c++;}}
static inline int hasedge(u64 code,int i,int j){return (code>>PIDX[i][j])&1ULL;}
static void code_to_adj(u64 code,int n,int*adj){
    for(int i=0;i<n;i++)adj[i]=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(hasedge(code,i,j)){adj[i]|=1<<j;adj[j]|=1<<i;}
}
static u64 adj_to_code(int n,const int*adj){
    u64 c=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if((adj[i]>>j)&1)c|=1ULL<<PIDX[i][j];return c;
}
static inline u64 mix(u64 h,u64 v){h^=v+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}
static int cmp_u64(const void*a,const void*b){u64 x=*(u64*)a,y=*(u64*)b;return x<y?-1:(x>y);}

static u64 invariant(const int *adj,int n){
    long long A[MAXN][MAXN],Pk[5][MAXN][MAXN];
    for(int i=0;i<n;i++)for(int j=0;j<n;j++) A[i][j]=(adj[i]>>j)&1;
    memcpy(Pk[0],A,sizeof A);
    for(int k=1;k<5;k++)
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){
            long long s=0; for(int t=0;t<n;t++) s+=Pk[k-1][i][t]*A[t][j];
            Pk[k][i][j]=s;
        }
    int col[MAXN]; for(int i=0;i<n;i++) col[i]=__builtin_popcount(adj[i]);
    for(int it=0;it<n;it++){
        u64 sig[MAXN];
        for(int v=0;v<n;v++){
            u64 nb[MAXN]; int c=0;
            for(int u=0;u<n;u++) if((adj[v]>>u)&1) nb[c++]=col[u];
            qsort(nb,c,8,cmp_u64);
            u64 h=1469598103934665603ULL; h=mix(h,col[v]);
            for(int i=0;i<c;i++) h=mix(h,nb[i]);
            sig[v]=h;
        }
        u64 s2[MAXN]; memcpy(s2,sig,8*n); qsort(s2,n,8,cmp_u64);
        int nc[MAXN],same=1;
        for(int v=0;v<n;v++) for(int k=0;k<n;k++) if(s2[k]==sig[v]){nc[v]=k;break;}
        for(int v=0;v<n;v++) if(nc[v]!=col[v]) same=0;
        memcpy(col,nc,sizeof(int)*n);
        if(same) break;
    }
    u64 vs[MAXN];
    for(int v=0;v<n;v++){
        u64 h=14695981039346656037ULL;
        h=mix(h,col[v]); h=mix(h,__builtin_popcount(adj[v]));
        for(int k=1;k<5;k++) h=mix(h,(u64)Pk[k][v][v]);
        u64 nb[MAXN]; int c=0;
        for(int u=0;u<n;u++) if((adj[v]>>u)&1) nb[c++]=(u64)(col[u]*1000+__builtin_popcount(adj[u]));
        qsort(nb,c,8,cmp_u64);
        for(int i=0;i<c;i++) h=mix(h,nb[i]);
        vs[v]=h;
    }
    u64 vsorted[MAXN]; memcpy(vsorted,vs,8*n); qsort(vsorted,n,8,cmp_u64);
    u64 ps[MAXN*MAXN]; int np=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
        u64 h=1099511628211ULL;
        u64 a=vs[i]<vs[j]?vs[i]:vs[j], b=vs[i]<vs[j]?vs[j]:vs[i];
        h=mix(h,a); h=mix(h,b);
        for(int k=0;k<5;k++) h=mix(h,(u64)Pk[k][i][j]);
        ps[np++]=h;
    }
    qsort(ps,np,8,cmp_u64);
    u64 H=1469598103934665603ULL; H=mix(H,n);
    for(int i=0;i<n;i++) H=mix(H,vsorted[i]);
    for(int i=0;i<np;i++) H=mix(H,ps[i]);
    return H;
}

typedef struct{u64*k;size_t cap,cnt;}HS;
static void hs_init(HS*h,size_t cap){h->cap=cap;h->cnt=0;h->k=calloc(cap,8);
    if(!h->k){fprintf(stderr,"out of memory (hash %zu)\n",cap);exit(1);} }
static int hs_add(HS*h,u64 key){
    u64 k=key|1; size_t i=(size_t)((k*11400714819323198485ULL)&(h->cap-1));
    while(h->k[i]){ if(h->k[i]==k) return 0; i=(i+1)&(h->cap-1); }
    h->k[i]=k; h->cnt++; return 1;
}
static int connected(const int*adj,int n){
    int seen=1,fr=1;
    while(fr){int nf=0;for(int v=0;v<n;v++)if((fr>>v)&1)nf|=adj[v];nf&=~seen;seen|=nf;fr=nf;}
    return seen==(1<<n)-1;
}
static void print_g6(const int*adj,int n){
    char out[64]; int o=0;
    out[o++]=(char)(n+63);
    int nb=n*(n-1)/2, acc=0, cnt=0;
    for(int j=1;j<n;j++) for(int i=0;i<j;i++){
        acc=(acc<<1)|((adj[i]>>j)&1);
        if(++cnt==6){ out[o++]=(char)(acc+63); acc=0; cnt=0; }
    }
    if(cnt){ acc<<=(6-cnt); out[o++]=(char)(acc+63); }
    (void)nb; out[o]=0;
    puts(out);
}

int main(int argc,char**argv){
    int nmax = argc>1 ? atoi(argv[1]) : 9;
    if(nmax<1||nmax>MAXN){fprintf(stderr,"n must be 1..%d\n",MAXN);return 1;}
    u64 *cur=malloc(8); cur[0]=0; size_t curN=1;
    for(int n=2;n<=nmax;n++){
        init_pidx(n-1);
        int (*adjs)[MAXN] = malloc(sizeof(int)*MAXN*curN);
        for(size_t g=0;g<curN;g++) code_to_adj(cur[g],n-1,adjs[g]);
        init_pidx(n);
        size_t want = (size_t)(A000088[n]*(u64)5/2);
        size_t cap=1; while(cap<want) cap<<=1;
        HS h; hs_init(&h,cap);
        u64 *out=malloc(8*(size_t)(A000088[n]+16));
        if(!out){fprintf(stderr,"out of memory (out)\n");return 1;}
        size_t on=0, nc=0;
        for(size_t g=0;g<curN;g++){
            for(int S=0;S<(1<<(n-1));S++){
                int adj[MAXN];
                for(int i=0;i<n-1;i++) adj[i]=adjs[g][i];
                adj[n-1]=S;
                for(int i=0;i<n-1;i++) if((S>>i)&1) adj[i]|=1<<(n-1);
                if(!hs_add(&h,invariant(adj,n))) continue;
                out[on++]=adj_to_code(n,adj);
                if(n==nmax && connected(adj,n)){ print_g6(adj,n); nc++; }
                else if(connected(adj,n)) nc++;
            }
        }
        free(adjs); free(h.k); free(cur);
        cur=out; curN=on;
        fprintf(stderr,"n=%2d: %9zu graphs (expected %9llu), %9zu connected (expected %9llu)%s\n",
                n,on,(unsigned long long)A000088[n],nc,(unsigned long long)A001349[n],
                (on==A000088[n]&&nc==A001349[n])?"  OK":"  *** MISMATCH ***");
        if(on!=A000088[n]||nc!=A001349[n]){
            fprintf(stderr,"de-duplication invariant was not complete; aborting.\n");
            return 2;
        }
    }
    return 0;
}
