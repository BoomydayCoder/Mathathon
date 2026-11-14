#!/usr/bin/env python3
"""
catalan_variants.py -- variants of the Hankel family for Catalan's constant G (functional II, 'cosh') and for
L(2,chi_-3) ('chi3'); EXACT finite-K computations (evidence only).  Columns as in catalan_table.py, per K^2.

 (A) nodes of one parity only (n = 2,4,..,2K or n = 1,3,..,2K-1): then the X-coefficients (-1)^n have constant
     sign and a reflection lemma holds (odd character, nodes = multiples of the conductor), the kernel is
     ~e^{-2 pi u} per unit node spacing (decay -> (9/2) log 2), but primes up to 4K enter.
 (B) Fauzan-type families  [mu_X(W t^{i+i'}/D_K)]_{i,i'<h},  W = prod_{n<=N}(t+n^2)^a, at K = 40
     (including the shape (K,N,a,h) = (40,3,6,37) of the zeta(5) proof), for 'cosh' and 'chi3'.
 (C) h = lam*K with W = 1.
"""
import sys, time
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from catalan_lib import *
import catalan_scan as CS
from multiprocessing import Pool

if __name__ == '__main__':
    print("(A) functional II on nodes of one parity")
    for (lab, j0) in (('even nodes n=2,4,..,2K', 1), ('odd nodes n=1,3,..,2K-1', 0)):
        for K in (20, 40, 60, 80):
            r = analyse(2, 'cosh', fmpq(1), K, j0=j0, step=2, want_vals=True)
            print(f"  {lab:26s} K={K:3d} deg={r['deg']:3d} sign={r['sign']:+d} cost={r['cost'] / K ** 2:.4f} "
                  f"logdet={r['logdet'] / K ** 2:.4f} total={r['total'] / K ** 2:+.4f}  v2={r['vals'].get(2, 0)}")
            if K == 60:
                print("      -v_p:", {p: -v for p, v in sorted(r['vals'].items())})
            sys.stdout.flush()
    K = 40
    cfgs = [(0, 0, 40), (0, 0, 37), (3, 6, 37), (3, 4, 37), (3, 2, 37), (2, 6, 38), (1, 6, 40), (4, 2, 36), (1, 2, 40), (2, 2, 40)]
    jobs = [(2, ker, K, 0, N, a, h) for ker in ('cosh', 'chi3') for (N, a, h) in cfgs]
    lams = (0.4, 0.55, 0.7, 0.8, 0.9, 1.0, 1.1, 1.3)
    jobs2 = [(2, 'cosh', KK, 0, 0, 0, int(round(l * KK))) for KK in (40, 60) for l in lams]
    with Pool(4) as pool:
        out = pool.map(CS.one, jobs, chunksize=1)
        out2 = pool.map(CS.one, jobs2, chunksize=1)
    print("(B) Fauzan-type families at K=40:  kernel N a h | deg cost logdet total   (per K^2)")
    for args, r in out:
        s, ker, K, j0, N, a, h = args
        print(f"  {ker:5s} N={N} a={a} h={h} | " + ("None" if r is None else f"deg={r[0]} cost={r[1]:.4f} logdet={r[2]:.4f} total={r[3]:+.4f}"))
    print("(C) functional II, W=1, h = lam*K:  K h | deg cost logdet total   (per K^2)")
    for args, r in out2:
        print(f"  K={args[2]:3d} h={args[6]:3d} | deg={r[0]:3d} cost={r[1]:.4f} logdet={r[2]:.4f} total={r[3]:+.4f}")
