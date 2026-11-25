// Builder B: Gram-based integer LLL (exact int64 Gram, long double GSO) and Fincke-Pohst enumeration.
#include <math.h>
#include <string.h>
#include <stdlib.h>
typedef long long ll;
#define GG(i,j) G[(size_t)(i)*n+(j)]
#define MU(i,j) mu[(size_t)(i)*n+(j)]
#define RR(i,j) r[(size_t)(i)*n+(j)]
long lll_gram(int n, ll *G, ll *U, double delta) {
    long double *mu = calloc((size_t)n*n, sizeof(long double));
    long double *r = calloc((size_t)n*n, sizeof(long double));
    int k = 1; long swaps = 0;
    RR(0,0) = (long double)GG(0,0);
    while (k < n) {
        for (int pass = 0; pass < 1000; pass++) {
            for (int j = 0; j <= k; j++) {
                long double s = (long double)GG(k,j);
                for (int i = 0; i < j; i++) s -= MU(j,i) * RR(k,i);
                RR(k,j) = s;
                if (j < k) MU(k,j) = s / RR(j,j);
            }
            int changed = 0;
            for (int j = k-1; j >= 0; j--) {
                long double q = roundl(MU(k,j));
                if (q != 0) {
                    ll qi = (ll)q; changed = 1;
                    for (int i = 0; i < n; i++) GG(k,i) -= qi*GG(j,i);
                    for (int i = 0; i < n; i++) GG(i,k) -= qi*GG(i,j);
                    for (int i = 0; i < n; i++) U[(size_t)k*n+i] -= qi*U[(size_t)j*n+i];
                    for (int i = 0; i < j; i++) MU(k,i) -= q*MU(j,i);
                    MU(k,j) -= q;
                }
            }
            if (!changed) break;
        }
        if (RR(k,k) < (delta - MU(k,k-1)*MU(k,k-1)) * RR(k-1,k-1)) {
            for (int i = 0; i < n; i++) { ll t = GG(k,i); GG(k,i) = GG(k-1,i); GG(k-1,i) = t; }
            for (int i = 0; i < n; i++) { ll t = GG(i,k); GG(i,k) = GG(i,k-1); GG(i,k-1) = t; }
            for (int i = 0; i < n; i++) { ll t = U[(size_t)k*n+i]; U[(size_t)k*n+i] = U[(size_t)(k-1)*n+i]; U[(size_t)(k-1)*n+i] = t; }
            swaps++;
            if (k > 1) k--; else { RR(0,0) = (long double)GG(0,0); k = 1; }
        } else k++;
    }
    free(mu); free(r);
    return swaps;
}

typedef struct { int n; const double *L; double bound; long maxc, cnt; int *out; double *vals; int x[200]; long nodes, maxnodes; int overflow; } ctx_t;
static void rec(ctx_t *c, int j, double partial, int allzero) {
    int n = c->n;
    if (c->cnt >= c->maxc || c->overflow) return;
    double s = 0.0;
    for (int i = j + 1; i < n; i++) s += c->x[i] * c->L[(size_t)i * n + j];
    double d = c->L[(size_t)j * n + j];
    double center = -s / d;
    double rem = c->bound - partial;
    if (rem < 0) return;
    double rad = sqrt(rem) / d;
    long lo = (long)ceil(center - rad - 1e-9), hi = (long)floor(center + rad + 1e-9);
    if (allzero && lo < 0) lo = 0;
    for (long v = lo; v <= hi; v++) {
        double t = v * d + s;
        double p = partial + t * t;
        if (p > c->bound * (1 + 1e-12)) continue;
        c->nodes++;
        if (c->nodes > c->maxnodes) { c->overflow = 1; return; }
        c->x[j] = (int)v;
        int az = allzero && (v == 0);
        if (j == 0) {
            if (az) continue;
            if (c->cnt < c->maxc) { memcpy(c->out + c->cnt * n, c->x, n * sizeof(int)); c->vals[c->cnt] = p; c->cnt++; }
            else return;
        } else {
            rec(c, j - 1, p, az);
            if (c->cnt >= c->maxc || c->overflow) return;
        }
    }
    c->x[j] = 0;
}
// enumerate all nonzero x (one of +-x) with |x L|^2 <= bound, where Gram = L L^T (L lower triangular)
long fp_enum(int n, const double *L, double bound, long maxc, int *out, double *vals, long maxnodes, int *flag, long *nodes) {
    ctx_t c; c.n = n; c.L = L; c.bound = bound; c.maxc = maxc; c.cnt = 0; c.out = out; c.vals = vals;
    c.nodes = 0; c.maxnodes = maxnodes; c.overflow = 0; memset(c.x, 0, sizeof(c.x));
    rec(&c, n - 1, 0.0, 1);
    *flag = c.overflow ? 1 : (c.cnt >= maxc ? 2 : 0);
    *nodes = c.nodes;
    return c.cnt;
}
