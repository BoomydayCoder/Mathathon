#!/usr/bin/env python3
"""
zeta5margin_exact_s.py -- exact valuations of Delta = det[mu_X^{(s)}(D_N^a t^{i+j}/D_K)]_{i,j<K-N} versus the
local cost MODEL (water filling in the class model with base a l_N - l_K - (s-1); generalised outer bound).
For s=5, a=6, N=3K/40 the generalised outer bound is asserted to coincide with Fauzan's (4.14).
Nothing is claimed by this file; it only computes.      Usage: python3 zeta5margin_exact_s.py s K N a
"""
import sys, time, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from hankel_lib import primes_upto, moments_F1, hankel, det_poly_affine, gauss_val, eval_poly_at_zeta, log_abs_arb
from zeta5margin_exact import gamma_out as gamma_out_415


def ell(A, c, p):
    return (A - c) // p + 1 + (A + c) // p


def gamma_wf(s, a, K, N, h, p):
    w = []
    for c in range(1, (p - 1) // 2 + 1):
        B = a * ell(N, c, p) - ell(K, c, p) - (s - 1)
        w += [B + 2 * i for i in range(h)]
    B0 = 2 * a * (N // p) - 2 * (K // p) + 1
    w += [B0 + 4 * i for i in range(h)]
    w.sort()
    return sum(w[:h])


def gamma_out(s, a, K, N, h, p):
    """generalisation of (4.11)-(4.14) to (s, a); needs p <= K < 3p, 2N < p, h = K - N."""
    if not (p >= s + 2 and p <= K < 3 * p and p * p > 2 * K and 2 * N < p and h == K - N):
        return None
    cost = 0; zeros = 0          # cost = -2 sum w (integer), zeros = number of zero weights
    for c in range(1, (p - 1) // 2 + 1):
        l, d = ell(K, c, p), ell(N, c, p)
        if d == 0:
            cost += (s + 2) * max(l - 2, 0); zeros += min(l, 2)
        else:
            ws2 = [min(0, 2 * i + a - (l + s - 1)) for i in range(max(l - 2, 0))]   # doubled weights
            cost += -sum(ws2); zeros += sum(1 for v in ws2 if v == 0) + (1 if l >= 2 else 0)
    mK = K // p
    cost += {1: 1, 2: 4}[mK]; zeros += {1: 0, 2: 1}[mK]
    r = max(0, 2 * h - K + a * N - (s - 1) * (p - 1) // 2)
    return -cost - min(r, zeros)


if __name__ == "__main__":
    s, K, N, a = (int(v) for v in sys.argv[1:5])
    h = K - N
    t0 = time.time()
    al, be = moments_F1(s, K, N, a, h)
    cs = det_poly_affine(hankel(al, h), hankel(be, h))
    val = eval_poly_at_zeta(cs, s, 4000 + 60 * K * K)
    logdet = log_abs_arb(val)[0]
    print(f"# s={s} K={K} N={N} a={a} h={h}: det computed in {time.time()-t0:.0f}s; log|Delta(zeta({s}))|/K^2 = {logdet/K**2:.5f}")
    print("#    p   K/p  | true |  out |   wf | model=max | true-model")
    tot_true = tot_model = 0.0
    viol = 0
    for p in primes_upto(2 * K + 20):
        tr = gauss_val(cs, p)
        valid = p >= s and p * p > 2.6 * K + 6 + 2 * a * N       # crude sufficient condition for the class model
        go = gamma_out(s, a, K, N, h, p)
        gw = gamma_wf(s, a, K, N, h, p) if valid else None
        if s == 5 and a == 6 and 40 * N == 3 * K and go is not None:
            assert go == gamma_out_415(K, N, p), (p, go, gamma_out_415(K, N, p))
        if p > K: model = 0
        else:
            cand = [v for v in (go, gw) if v is not None]
            model = max(cand) if cand else None
        tot_true += -tr * math.log(p)
        if model is not None:
            tot_model += -model * math.log(p)
            if model > tr: viol += 1
        else:
            tot_model += -tr * math.log(p)          # small primes: no model, use the truth
        if p <= K or tr != 0:
            f = lambda v: '    .' if v is None else f"{v:5d}"
            print(f" {p:5d} {K/p:6.3f} | {tr:5d}|{f(go)} |{f(gw)} | {f(model)}     | {f(None if model is None else tr - model)}")
    print(f"# model > true at {viol} primes")
    print(f"# K^-2 * [ log|Delta(zeta)| , cost_true = -sum v_p log p , cost_model (truth at small primes) ] = {logdet/K**2:.5f}, {tot_true/K**2:.5f}, {tot_model/K**2:.5f}")
    print(f"# finite-K totals: true {(logdet + tot_true)/K**2:+.5f}   model {(logdet + tot_model)/K**2:+.5f}")
    print(f"# time {time.time()-t0:.0f}s")
