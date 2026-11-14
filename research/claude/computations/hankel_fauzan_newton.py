#!/usr/bin/env python3
"""
hankel_fauzan_newton.py -- does the "two-sided Newton basis" also work for Fauzan's zeta(5) family?

Family: K=40n, N=3n, h=37n=K-N, W = D_N^6, s=5.  Bases (unitriangular, so det is unchanged):
    rows  f_i  = prod_{j=N+1}^{N+i} (t+j^2),   columns g_i' = prod_{j=K-i'+1}^{K} (t+j^2).
Entries: mu_X( D_N^5 / prod_{j=N+i+1}^{K-i'}(t+j^2) )  resp.  mu( D_N^5 prod_{j=K-i'+1}^{N+i}(t+j^2) ).
We compare the true Gauss valuation v_p(Delta_K) with the tropical determinant (min-cost perfect matching)
of the matrix of entrywise valuations, which is the best bound obtainable from row/column weights.
Finding (n=1,2): equality for the large primes (p >= 23 at K=40, p >= 47 at K=80), a loss for smaller p
(total loss about 0.5 K^2), i.e. this basis does NOT capture the gain from the zeros of W.
Usage: python3 hankel_fauzan_newton.py n [s]
"""
import sys, time, math, pickle
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import *
import numpy as np
from scipy.optimize import linear_sum_assignment

def mu_rational(num, lo, hi, s, H, mu_mon):
    """mu_X( num(t) / prod_{j=lo}^{hi} (t+j^2) ) -> (alpha, beta); num fmpz_poly. If lo>hi, denominator = 1."""
    if lo > hi:
        cs = num.coeffs()
        return sum((fmpq(int(cs[k]))*mu_mon[k] for k in range(len(cs))), fmpq(0)), fmpq(0)
    den = fmpz_poly([1])
    for j in range(lo, hi+1):
        den = den*fmpz_poly([j*j, 1])
    al = fmpq(0); be = fmpq(0)
    if num.degree() >= den.degree():
        quo, rem = divmod(num, den)
        cs = quo.coeffs()
        al += sum((fmpq(int(cs[k]))*mu_mon[k] for k in range(len(cs))), fmpq(0))
    dp = den.derivative()
    for j in range(lo, hi+1):
        nj = num(fmpz(-j*j))
        if nj == 0: continue
        r = fmpq(nj, dp(fmpz(-j*j)))
        be += r*(j**(s-1))
        al += r*(-(j**(s-1))*H[j] - fmpq(1, s-1) + fmpq(1, 2*j))
    return al, be

def newton_fauzan(s, K, N, a):
    h = K - N
    H = harmonic_s(K, s)
    mu_mon = [mu_monomial(e, s) for e in range(0, (a-1)*N + h + 5)]
    DN = D_poly(N)
    W = DN**(a-1)
    Al = [[None]*h for _ in range(h)]; Be = [[None]*h for _ in range(h)]
    for i in range(h):
        for ip in range(h):
            lo, hi = N+i+1, K-ip
            if lo <= hi:
                Al[i][ip], Be[i][ip] = mu_rational(W, lo, hi, s, H, mu_mon)
            else:
                P = W
                for j in range(hi+1, lo):
                    P = P*fmpz_poly([j*j, 1])
                Al[i][ip], Be[i][ip] = mu_rational(P, 1, 0, s, H, mu_mon)
    return Al, Be

def ve(a, b, p):
    v = 10**5
    if a != 0: v = min(v, vp(a, p))
    if b != 0: v = min(v, vp(b, p))
    return v

if __name__ == "__main__":
    n = int(sys.argv[1]); s = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    K, N, a = 40*n, 3*n, 6
    h = K - N
    t0 = time.time()
    Al, Be = newton_fauzan(s, K, N, a)
    A = fmpq_mat(h, h, [Al[i][j] for i in range(h) for j in range(h)])
    B = fmpq_mat(h, h, [Be[i][j] for i in range(h) for j in range(h)])
    cs = det_poly_affine(A, B)
    al, be = moments_F1(s, K, N, a, h)
    cs0 = det_poly_affine(hankel(al, h), hankel(be, h))
    print(f"s={s} K={K} N={N} h={h}: det(Newton tail basis) == det(monomial): {all(x == y for x, y in zip(cs, cs0))}  ({time.time()-t0:.0f}s)")
    tot_true = 0.0; tot_trop = 0.0
    print("   p | true v_p(Delta_K) | tropical det (Newton tail basis)")
    for p in primes_upto(2*K):
        V = np.zeros((h, h))
        for i in range(h):
            for j in range(h):
                V[i, j] = ve(Al[i][j], Be[i][j], p)
        r, c = linear_sum_assignment(V)
        trop = int(V[r, c].sum()); true = gauss_val(cs, p)
        tot_true += true*math.log(p); tot_trop += trop*math.log(p)
        if true != 0 or trop != 0:
            print(f"   {p:3d} | {true:6d} | {trop:6d}")
    print(f"   sum v_p log p / K^2:  true {tot_true/K**2:.4f}   tropical bound {tot_trop/K**2:.4f}   (difference {(tot_true-tot_trop)/K**2:.4f})")
