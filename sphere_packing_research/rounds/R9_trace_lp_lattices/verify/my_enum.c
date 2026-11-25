// Independent enumeration (verifier's own code): all nonzero x in Z^n, one per +-pair,
// with |R x|^2 <= bound, R upper triangular (row-major n*n, G = R^T R).
#include <math.h>
#include <string.h>
static int n_; static const double *R_; static double bound_;
static long maxc_, cnt_, nodes_, maxnodes_; static int *out_; static double *vals_;
static int x_[64]; static int stop_;
static void go(int k, double acc, int lead) {
    double s = 0.0;
    for (int j = k + 1; j < n_; j++) s += R_[k * n_ + j] * x_[j];
    double d = R_[k * n_ + k];
    double rem = bound_ - acc;
    if (rem < 0) return;
    double rad = sqrt(rem) / fabs(d);
    double c = -s / d;
    long lo = (long)ceil(c - rad - 1e-7), hi = (long)floor(c + rad + 1e-7);
    if (lead && lo < 0) lo = 0;
    for (long v = lo; v <= hi; v++) {
        if (stop_) return;
        double y = d * (double)v + s;
        double a2 = acc + y * y;
        if (a2 > bound_) continue;
        if (++nodes_ > maxnodes_) { stop_ = 1; return; }
        x_[k] = (int)v;
        int l2 = lead && (v == 0);
        if (k == 0) {
            if (l2) continue;
            if (cnt_ >= maxc_) { stop_ = 2; return; }
            memcpy(out_ + cnt_ * n_, x_, n_ * sizeof(int));
            vals_[cnt_] = a2; cnt_++;
        } else go(k - 1, a2, l2);
    }
    x_[k] = 0;
}
long enum_short(int n, const double *R, double bound, long maxc, int *out, double *vals,
                long maxnodes, int *status, long *nodes) {
    n_ = n; R_ = R; bound_ = bound; maxc_ = maxc; cnt_ = 0; nodes_ = 0; maxnodes_ = maxnodes;
    out_ = out; vals_ = vals; stop_ = 0; memset(x_, 0, sizeof(x_));
    go(n - 1, 0.0, 1);
    *status = stop_; *nodes = nodes_;
    return cnt_;
}
