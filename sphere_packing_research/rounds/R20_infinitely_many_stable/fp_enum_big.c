// Fincke-Pohst enumeration of all nonzero integer vectors x (one of each +-pair)
// with x^T G x <= bound, where G = L L^T (L lower triangular, row-major n*n).
// x^T G x = sum_j ( sum_{i>=j} x_i L[i][j] )^2
#include <math.h>
#include <string.h>

typedef struct {
    int n;
    const double *L;
    double bound;
    long maxc;
    long cnt;
    int *out;
    double *vals;
    int x[1024];
    long nodes;
    long maxnodes;
    int overflow;
} ctx_t;

static void rec(ctx_t *c, int j, double partial, int allzero) {
    int n = c->n;
    if (c->cnt >= c->maxc || c->overflow) return;
    double s = 0.0;
    for (int i = j + 1; i < n; i++) s += c->x[i] * c->L[i * n + j];
    double d = c->L[j * n + j];
    double center = -s / d;
    double rem = c->bound - partial;
    if (rem < 0) return;
    double r = sqrt(rem) / d;
    long lo = (long)ceil(center - r - 1e-9), hi = (long)floor(center + r + 1e-9);
    if (allzero && lo < 0) lo = 0;
    for (long v = lo; v <= hi; v++) {
        double t = v * d + s;
        double p = partial + t * t;
        if (p > c->bound) continue;
        c->nodes++;
        if (c->nodes > c->maxnodes) { c->overflow = 1; return; }
        c->x[j] = (int)v;
        int az = allzero && (v == 0);
        if (j == 0) {
            if (az) continue;
            if (c->cnt < c->maxc) {
                memcpy(c->out + c->cnt * n, c->x, n * sizeof(int));
                c->vals[c->cnt] = p;
                c->cnt++;
            } else return;
        } else {
            rec(c, j - 1, p, az);
            if (c->cnt >= c->maxc || c->overflow) return;
        }
    }
}

// returns count (>=0); sets *flag = 1 if node limit hit, 2 if maxc reached
long fp_enum(int n, const double *L, double bound, long maxc, int *out, double *vals,
             long maxnodes, int *flag, long *nodes) {
    ctx_t c;
    c.n = n; c.L = L; c.bound = bound; c.maxc = maxc; c.cnt = 0; c.out = out; c.vals = vals;
    c.nodes = 0; c.maxnodes = maxnodes; c.overflow = 0;
    memset(c.x, 0, sizeof(c.x));
    rec(&c, n - 1, 0.0, 1);
    *flag = c.overflow ? 1 : (c.cnt >= maxc ? 2 : 0);
    *nodes = c.nodes;
    return c.cnt;
}
