#!/usr/bin/env python3
"""
hankel_theoremA_check.py -- end-to-end check of the arithmetic theorem for s=3, W=1, h=K:

   Phi_K * det G_K(X) in Z[X],   Phi_K = prod_p p^{e_p},
   e_p = 7K floor(log_p 3K)          for p <= P0 := sqrt(2K)+3      (crude bound, all small primes)
       = 3K - 2 floor(K/p)           for P0 < p <= K/2
       = 6(K-p+1)                    for max(P0, K/2) < p <= K
       = 0                           for p > K.
For each K we verify v_p(det G_K) >= -e_p for every prime (Gauss valuation, exact arithmetic) and that
no other prime occurs in a denominator.  We also print, on the scale K^2:
   true cost = -log content(det G_K);   main term M_K = [sum_{p<=K/2} 3K log p + sum_{K/2<p<=K} 6(K-p) log p]/K^2
   (-> 9/4 by the prime number theorem; = 6 int_{K/2}^{K} theta(u) du / K^2);
   log Phi_K / K^2 (dominated at these K by the crude small-prime factor, which is O(K^{3/2}));
   log det G_K(zeta(3)) / K^2  (-> -(9/2) log 2 = -3.1192 conjecturally; proved: <= -(1+2 log(5/2)) = -2.8326 + o(1)).
Extra column (NOT proved for the primes p <= sqrt(2K)+3): is  lcm(1..K)^{3K} * det G_K  in Z[X] ?
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import *

def exponents(K):
    P0 = math.sqrt(2 * K) + 3
    e = {}
    for p in primes_upto(K):
        if p <= P0:
            e[p] = 7 * K * int(math.floor(math.log(3 * K) / math.log(p) + 1e-12))
        elif 2 * p <= K:
            e[p] = 3 * K - 2 * (K // p)
        else:
            e[p] = 6 * (K - p + 1)
    return e

if __name__ == '__main__':
    Ks = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [10, 20, 30, 41, 60, 80, 100]
    print("  K | Phi_K det in Z[X] | other primes in denominators | true cost/K^2 | M_K | log Phi_K/K^2 | log det(zeta3)/K^2 | log|prim(zeta3)|/K^2 | d_K^{3K} det in Z[X]")
    for K in Ks:
        al, be = moments_F1(3, K, 0, 0, K)
        cs = det_poly_affine(hankel(al, K), hankel(be, K))
        e = exponents(K)
        ok = True
        ok_lcm = True
        cont = content(cs)
        rest = fmpq(1, cont.q)   # denominator part of the content
        for p in primes_upto(max(3 * K, 50)):
            v = gauss_val(cs, p)
            if v < -e.get(p, 0):
                ok = False
            elcm = 3 * K * int(math.floor(math.log(K) / math.log(p) + 1e-12)) if p <= K else 0
            if v < -elcm:
                ok_lcm = False
        den = int(cont.q)
        for p in primes_upto(K):
            while den % p == 0:
                den //= p
        cost = -log_abs_fmpq(cont)
        M = sum(min(3 * K, 6 * (K - p)) * math.log(p) for p in primes_upto(K)) / K ** 2
        logPhi = sum(ep * math.log(p) for p, ep in e.items()) / K ** 2
        val = eval_poly_at_zeta(cs, 3, 6000 + 80 * K * K)
        ld, _ = log_abs_arb(val)
        print(f"  {K:3d} | {ok} | {'none' if den == 1 else den} | {cost / K ** 2:.4f} | {M:.4f} | {logPhi:.4f} | {ld / K ** 2:.4f} | {(ld + cost) / K ** 2:+.4f} | {ok_lcm and den == 1}")
        sys.stdout.flush()
