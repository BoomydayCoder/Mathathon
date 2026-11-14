#!/usr/bin/env python3
"""
zeta5margin_exact.py -- exact check of the water-filling (class-model) bound for Fauzan's zeta(5)
determinant at ALL primes, in particular in the outer range p > K/3.

  Delta_K = det[mu_X(D_N^6 t^{i+j}/D_K)]_{i,j<h},  K=40n, N=3n, h=37n   (Fauzan (2.4)).
  true_p   = exact Gauss valuation v_p^G(Delta_K)                        (python-flint, exact)
  out_p    = Fauzan's gamma_p^out, formula (4.14), when (4.9) holds
  wf_p     = 2 * (sum of the h smallest elements of the multiset
               { (6 l_N(c) - l_K(c) - 4)/2 + i : 1<=c<=(p-1)/2, i>=0 }  u
               { (12 m_N - 2 m_K + 1)/2 + 2 i : i>=0 } )
             = optimum of the class model (4.2)-(4.3) over allocations ("water filling")
  new_p    = max(out_p, wf_p)
Claim tested: wf_p <= true_p for every prime p >= 5 with p^2 > 2.6 K + 6   (Theorem W of the report).
Note: for p > K ('beyond') the applicable bound is 0 (Prop. 4.3); the wf column there is informational only
(still a valid lower bound, but useless).  'small' = primes where Theorem W is not applicable.
Nothing is claimed by this file; it only computes.     Usage: python3 zeta5margin_exact.py n
"""
import sys, time, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import primes_upto, moments_F1, hankel, det_poly_affine, gauss_val


def ell(A, c, p):
    return (A - c) // p + 1 + (A + c) // p


def gamma_wf(K, N, h, p):
    """doubled weights: class c>=1: B_c + 2i, B_c = 6 l_N - l_K - 4; class 0: B_0 + 4i, B_0 = 12 m_N - 2 m_K + 1."""
    w = []
    for c in range(1, (p - 1) // 2 + 1):
        B = 6 * ell(N, c, p) - ell(K, c, p) - 4
        w += [B + 2 * i for i in range(h)]
    B0 = 12 * (N // p) - 2 * (K // p) + 1
    w += [B0 + 4 * i for i in range(h)]
    w.sort()
    return sum(w[:h]), w[h - 1]          # (gamma_wf, doubled top level)


def gamma_out(K, N, p):
    if not (p >= 7 and p <= K < 3 * p and p * p > 2 * K and 2 * N < p and 5 * N <= 2 * p - 2):
        return None
    v = K - p * (K // p)
    u = max(0, N + v - p + 1)
    tp = min(N, v) + u
    rp = max(0, K + 4 * N - 2 * p + 2)
    if K < 2 * p:
        return -7 * (K - p) + 6 * tp - 1 - min(rp, p - 1 - N + u)
    return -7 * (K - p) + 3 + 12 * N + 5 * tp - min(rp, p + u)


if __name__ == "__main__":
    n = int(sys.argv[1])
    K, N = 40 * n, 3 * n
    h = K - N
    t0 = time.time()
    al, be = moments_F1(5, K, N, 6, h)
    cs = det_poly_affine(hankel(al, h), hankel(be, h))
    print(f"# s=5 K={K} N={N} h={h}: Delta_K computed exactly in {time.time()-t0:.0f}s")
    P = primes_upto(K + 30)
    print("#    p   K/p    | true |  out |   wf | new=max | true-new | 2*level | wf valid (p>=5, p^2>2.6K+6) | wf<=true")
    bad = 0
    loss = {'out': 0.0, 'new': 0.0, 'wf': 0.0}
    lossr = {}
    for p in P:
        tr = gauss_val(cs, p)
        go = gamma_out(K, N, p)
        valid = p >= 5 and p * p > 104 * n + 6
        gw, lev = gamma_wf(K, N, h, p) if p >= 3 else (None, None)
        cand = [x for x in (go, gw if valid else None) if x is not None]
        new = max(cand) if cand else None
        ok = (gw <= tr) if valid else None
        if valid and not ok:
            bad += 1
        f = lambda x: '    .' if x is None else f"{x:5d}"
        rng = 'small' if not valid else ('inner' if 3 * p <= K else ('outer' if p <= K else 'beyond'))
        print(f" {p:5d} {K/p:6.3f} {rng:6s}| {tr:5d}|{f(go)} |{f(gw)} | {f(new)}   | {f(None if new is None else tr-new)}    | {f(lev)}   | {str(valid):5s} | {ok}")
        if new is not None:
            lossr.setdefault(rng, [0.0, 0.0, 0.0])
            lossr[rng][0] += (tr - new) * math.log(p) / K ** 2
            if go is not None:
                lossr[rng][1] += (tr - go) * math.log(p) / K ** 2
                lossr[rng][2] += (new - go) * math.log(p) / K ** 2
    print(f"# violations of wf_p <= true_p among valid primes: {bad}")
    for rng, (a, b, c) in lossr.items():
        print(f"# range {rng:6s}: sum (true-new) log p /K^2 = {a:.5f};  sum (true-out) log p/K^2 = {b:.5f};  gain sum (new-out) log p/K^2 = {c:.5f}")
    print(f"# total time {time.time()-t0:.0f}s")
