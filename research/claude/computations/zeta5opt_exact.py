#!/usr/bin/env python3
"""
zeta5opt_exact.py -- exact check of the GENERALISED valuation bounds for the family
    Delta = det[ mu_X( D_N^a t^{i+j} / D_K ) ]_{i,j<h}        (s = 5, a >= 1 integer, h <= K - N or any h).

Bounds tested against the exact Gauss valuation v_p^G(Delta) (python-flint, exact rationals):
  wf_p  (Theorem W_a): 2 * (sum of the h smallest elements of the multiset
           { (a l_N(c) - l_K(c) - 4)/2 + i : 1<=c<=(p-1)/2, i>=0 } u { (2a m_N - 2 m_K + 1)/2 + 2i : i>=0 } ),
        claimed for primes p >= 5 with p^2 > max(2K, 4h + 2aN - 2K + 2).           [any h]
  out_p (generalised Fauzan (4.12)-(4.14), computed class by class, NOT from a closed formula):
        classes c=1..(p-1)/2 with l = l_K(c), delta = [c <= N]:
            w_i = min(0, i + a*delta/2 - (l+4)/2)  (i < l-2),   w_i = 0  (l-2 <= i < l-delta);
        class 0: one pole -> w = -1/2;  two poles -> w = -2, 0;
        out_p = 2 sum w - min(r_p, #zero weights),  r_p = max(0, 2h + aN - K - 2p + 2),
        claimed for h = K - N and  p >= 7, p <= K < 3p, p^2 > 2K, 2N < p.
  zero_p: v_p >= 0 for p > K  (needs 2h + aN - K + 1 <= 2p; reported).
For a = 6 the closed formula (4.14) of Fauzan is compared with out_p (must agree).
Usage: python3 zeta5opt_exact.py a N K [h]
Nothing is claimed by this file; it only computes.
"""
import sys, time, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import primes_upto, moments_F1, hankel, det_poly_affine, gauss_val


def ell(A, c, p):
    return (A - c) // p + 1 + (A + c) // p


def gamma_wf(a, K, N, h, p):
    w = []
    for c in range(1, (p - 1) // 2 + 1):
        B = a * ell(N, c, p) - ell(K, c, p) - 4
        w += [B + 2 * i for i in range(h)]
    B0 = 2 * a * (N // p) - 2 * (K // p) + 1
    w += [B0 + 4 * i for i in range(h)]
    w.sort()
    return sum(w[:h])


def wf_valid(a, K, N, h, p):
    return p >= 5 and p * p > max(2 * K, 4 * h + 2 * a * N - 2 * K + 2)


def gamma_out(a, K, N, h, p):
    if h != K - N or not (p >= 7 and p <= K < 3 * p and p * p > 2 * K and 2 * N < p):
        return None
    tw = 0          # doubled weights
    z = 0
    for c in range(1, (p - 1) // 2 + 1):
        l = ell(K, c, p)
        d = 1 if c <= N else 0
        for i in range(l - d):
            if i < l - 2:
                w2 = min(0, 2 * i + a * d - (l + 4))
            else:
                w2 = 0
            tw += w2
            z += (w2 == 0)
    mK = K // p
    if mK == 1:
        tw += -1
    elif mK == 2:
        tw += -4
        z += 1
    r = max(0, 2 * h + a * N - K - 2 * p + 2)
    return tw - min(r, z)


def fauzan_414(K, N, p):
    v = K - p * (K // p)
    u = max(0, N + v - p + 1)
    tp = min(N, v) + u
    rp = max(0, K + 4 * N - 2 * p + 2)
    if K < 2 * p:
        return -7 * (K - p) + 6 * tp - 1 - min(rp, p - 1 - N + u)
    return -7 * (K - p) + 3 + 12 * N + 5 * tp - min(rp, p + u)


if __name__ == "__main__":
    a, N, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    h = int(sys.argv[4]) if len(sys.argv) > 4 else K - N
    t0 = time.time()
    al, be = moments_F1(5, K, N, a, h)
    cs = det_poly_affine(hankel(al, h), hankel(be, h))
    deg = max(i for i, c in enumerate(cs) if c != 0)
    print(f"# s=5 a={a} K={K} N={N} h={h}: Delta computed exactly in {time.time()-t0:.0f}s; deg_X = {deg}")
    P = primes_upto(max(K, 2 * h) + 30)
    print("#    p   K/p   | true |  wf  | out  | best | true-best | wf valid | wf ok | out ok")
    bad_wf = bad_out = bad_zero = n_wf = n_out = n_zero = mism = 0
    for p in P:
        tr = gauss_val(cs, p)
        valid = wf_valid(a, K, N, h, p)
        gw = gamma_wf(a, K, N, h, p) if p >= 3 else None
        go = gamma_out(a, K, N, h, p)
        if a == 6 and go is not None and 5 * N <= 2 * p - 2 and go != fauzan_414(K, N, p):
            mism += 1
        okw = oko = None
        cand = []
        if valid:
            n_wf += 1; okw = gw <= tr; bad_wf += (not okw); cand.append(gw)
        if go is not None:
            n_out += 1; oko = go <= tr; bad_out += (not oko); cand.append(go)
        if p > K:
            n_zero += 1; bad_zero += (tr < 0)
            if 2 * h + a * N - K + 1 <= 2 * p: cand.append(0)
        best = max(cand) if cand else None
        f = lambda x: '    .' if x is None else f"{x:5d}"
        print(f" {p:5d} {K/p:6.3f} | {tr:5d}|{f(gw)} |{f(go)} |{f(best)} | {f(None if best is None else tr-best)}     | {str(valid):5s}    | {okw} | {oko}")
    print(f"# Theorem W_a: {n_wf} primes tested, violations {bad_wf}")
    print(f"# generalised outer bound: {n_out} primes tested, violations {bad_out}")
    print(f"# v_p >= 0 for p > K: {n_zero} primes tested, violations {bad_zero}")
    if a == 6:
        print(f"# class-by-class out_p vs closed formula (4.14) (a=6): mismatches {mism}")
    print(f"# total time {time.time()-t0:.0f}s")
