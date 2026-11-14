#!/usr/bin/env python3
"""
zeta5margin_primesum.py -- finite-K sanity check of the new outer constant:
  S_old(K) = K^-2 sum_{K/3<p<=2h} (-v_p(S_K) - gamma_p^out) log p          -> I_out      = 127751/96000  = 1.330740
  S_new(K) = K^-2 sum_{K/3<p<=2h} (-v_p(S_K) - max(gamma^out,gamma^wf)) log p -> I_out^new = 2139737/1632000 = 1.311113
and of gamma_p^wf = p*Gamma_wf(K/p) + O(1) (discrete water filling vs. exact continuum).
Nothing is claimed by this file; it only computes.     Usage: python3 zeta5margin_primesum.py n [n ...]
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
import numpy as np
from fractions import Fraction as Fr
from zeta5margin_exact import gamma_out
from zeta5margin_asymptotic import Model, AL, LA

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def gamma_wf_np(K, N, h, p):
    c = np.arange(1, (p - 1) // 2 + 1)
    B = 6 * ((N - c) // p + 1 + (N + c) // p) - ((K - c) // p + 1 + (K + c) // p) - 4
    imax = (2 * h) // len(c) + 30
    i = np.arange(imax)
    w = (B[:, None] + 2 * i[None, :]).ravel()
    w = np.concatenate([w, 12 * (N // p) - 2 * (K // p) + 1 + 4 * i])
    w = np.partition(w, h - 1)[:h]
    assert w.max() < B.min() + 2 * (imax - 1)
    return int(w.sum())

def vS(K, N, h, p):   # Legendre (5.3), p odd, p^2 > 2h
    return 2 * h * (K // p) - 12 * h * (N // p) - 2 * sum((2 * i) // p for i in range(1, h))

if __name__ == "__main__":
    F = Model(5, 6, AL, LA)
    for n in (int(v) for v in sys.argv[1:]):
        K, N, h = 40 * n, 3 * n, 37 * n
        So = Sn = 0.0; maxdev = 0.0; worst = None; npr = 0
        for p in primes_upto(2 * h):
            if 3 * p <= K: continue
            if p > K: g_o = g_n = 0
            else:
                g_o = gamma_out(K, N, p); g_w = gamma_wf_np(K, N, h, p)
                g_n = max(g_o, g_w)
                if K < 3 * p and 2 * p < K and npr % 7 == 0:      # sample the range 2 < x < 3
                    dev = g_w - p * F.wf(Fr(K, p))[0]
                    if abs(dev) > maxdev: maxdev, worst = abs(float(dev)), p
                npr += 1
            v = vS(K, N, h, p)
            So += (-v - g_o) * math.log(p); Sn += (-v - g_n) * math.log(p)
        print(f"K={K}: S_old = {So/K**2:.6f} (limit 1.330740)   S_new = {Sn/K**2:.6f} (limit 1.311113)   gain = {(So-Sn)/K**2:.6f} (limit 0.019626);"
              f"   max |gamma_wf - p*Gamma_wf(K/p)| over sampled primes in (K/3,K/2): {maxdev:.3f} (p={worst})", flush=True)
