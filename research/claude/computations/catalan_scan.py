#!/usr/bin/env python3
"""
catalan_scan.py -- does a Fauzan-type modification rescue the Hankel method for Catalan's constant?

Family  F(K, j0, N, a, h):   G(X) = [ mu_X( W(t) t^{i+i'} / D(t) ) ]_{0<=i,i'<h},
   D(t) = prod_{n=j0+1}^{j0+K} (t + n^2),   W(t) = prod_{n=j0+1}^{j0+N} (t + n^2)^a,
for functional II ('cosh': integer nodes, weight (pi y/4) sinh(pi y)/cosh(pi y)^2, X0 = G)
and, for comparison, 'bose' s=3 (zeta(3)) and s=5.

Output per parameter set: degree in X, cost = -log(content), logdet = log|det(X0)|, total = log|primitive poly at X0|,
all divided by K^2.  total < 0 (growing like K^2) is what an irrationality proof needs.
EXACT computations for finite K: evidence only.
Usage: python3 catalan_scan.py full K [kernel]        full (N,a,h,j0) scan at this K; kernel in {cosh, chi3}
       python3 catalan_scan.py lam kernel K1 K2 ..    only W=1, j0=0, h = lam*K for lam in 0.4..1.3
"""
import sys, math, time
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from catalan_lib import *
from flint import fmpz_poly
from multiprocessing import Pool


def moments_W(s, kernel, K, j0, N, a, h):
    ns = list(range(j0 + 1, j0 + K + 1))
    pv_all = pole_values(s, kernel, 1, j0 + K)
    D = fmpz_poly([1])
    for n in ns:
        D = D * fmpz_poly([n * n, 1])
    W = fmpz_poly([1])
    for n in ns[:N]:
        W = W * fmpz_poly([n * n, 1]) ** a
    Dp = D.derivative()
    res = {}
    for n in ns:
        wv = W(fmpz(-n * n))
        if wv != 0:
            res[n] = fmpq(wv, Dp(fmpz(-n * n)))
    nmom = 2 * h - 1
    maxdeg = max(0, W.degree() + nmom - 1 - K)
    mm = [mu_monomial(e, s, kernel) for e in range(maxdeg + 1)]
    al, be = [], []
    Wt = W
    for e in range(nmom):
        a0 = fmpq(0)
        b0 = fmpq(0)
        if Wt.degree() >= K:
            quo, rem = divmod(Wt, D)
            cs = quo.coeffs()
            a0 += sum((fmpq(int(cs[k])) * mm[k] for k in range(len(cs))), fmpq(0))
        for n, r in res.items():
            rj = r * fmpq((-n * n) ** e)
            a0 += rj * pv_all[n - 1][0]
            b0 += rj * pv_all[n - 1][1]
        al.append(a0)
        be.append(b0)
        Wt = Wt * fmpz_poly([0, 1])
    return al, be


def X0(s, kernel, prec):
    return X0_arb(s, kernel, fmpq(1), prec)          # cosh: Catalan;  chi3: L(2,chi_-3)


def one(args):
    s, kernel, K, j0, N, a, h = args
    al, be = moments_W(s, kernel, K, j0, N, a, h)
    cs = det_poly_affine(hankel(al, h), hankel(be, h))
    nz = [k for k, c in enumerate(cs) if c != 0]
    if not nz:
        return args, None
    deg = max(nz)
    cont = content(cs)
    prec = 2000 + 20 * (K + j0 + a * N) ** 2
    while True:
        old = ctx.prec
        ctx.prec = prec
        val = eval_poly(cs, X0(s, kernel, prec))
        ctx.prec = old
        if (not val.contains(0)) and val.rad() < abs(val.mid()) * 1e-20:
            break
        prec *= 2
        if prec > 10 ** 7:
            return args, None
    ld, _ = log_abs_arb(val)
    cost = -log_abs_fmpq(cont)
    return args, (deg, cost / K ** 2, ld / K ** 2, (ld + cost) / K ** 2)


def lam_scan(kernel, Ks):
    lams = (0.4, 0.55, 0.7, 0.8, 0.9, 1.0, 1.1, 1.3)
    jobs = [(2, kernel, K, 0, 0, 0, int(round(l * K))) for K in Ks for l in lams]
    with Pool(4) as pool:
        out = pool.map(one, jobs, chunksize=1)
    print(f"# kernel {kernel}, W=1, j0=0, h = lam*K:  K h | deg cost/K^2 logdet/K^2 total/K^2")
    for args, r in out:
        print(f"  K={args[2]:3d} h={args[6]:3d} | deg={r[0]:3d} cost={r[1]:.4f} logdet={r[2]:.4f} total={r[3]:+.4f}")


def full_scan(kernel, K):
    jobs = []
    for (a, Nf) in [(0, 0)] + [(a, f) for a in (2, 4, 6, 10) for f in (0.05, 0.1, 0.2, 0.3)]:
        N = int(round(Nf * K))
        for lam in (0.7, 0.85, 1.0, 1.15, 1.3):
            jobs.append((2, kernel, K, 0, N, a, int(round(lam * K))))
    for j0f in (0.25, 0.5, 1.0):
        for lam in (0.85, 1.0, 1.15):
            jobs.append((2, kernel, K, int(round(j0f * K)), 0, 0, int(round(lam * K))))
    jobs = sorted(set(jobs))
    t = time.time()
    with Pool(4) as pool:
        out = pool.map(one, jobs, chunksize=1)
    print(f"# K={K}  kernel {kernel}; columns: j0 N a h | deg cost/K^2 logdet/K^2 total/K^2   ({time.time() - t:.0f}s)")
    best = None
    for args, r in out:
        s, kernel, K, j0, N, a, h = args
        if r is None:
            print(f"  j0={j0:3d} N={N:3d} a={a:2d} h={h:3d} | det = 0 or not evaluated")
            continue
        print(f"  j0={j0:3d} N={N:3d} a={a:2d} h={h:3d} | deg={r[0]:3d} cost={r[1]:.4f} logdet={r[2]:.4f} total={r[3]:+.4f}")
        if best is None or r[3] < best[1][3]:
            best = (args, r)
    print("# best total/K^2:", best)


if __name__ == '__main__':
    # usage: catalan_scan.py full K [kernel]   |   catalan_scan.py lam kernel K1 K2 ...
    if sys.argv[1] == 'lam':
        lam_scan(sys.argv[2], [int(x) for x in sys.argv[3:]])
    else:
        full_scan(sys.argv[3] if len(sys.argv) > 3 else 'cosh', int(sys.argv[2]))
