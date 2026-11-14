#!/usr/bin/env python3
"""
anand_zeta7_constants.py -- independent recomputation of the arithmetic constants of
P. Anand, "zeta(7) is irrational" (Zenodo 22920911), from the preprint's OWN definitions.

Parts
 A. Outer integral I_out (73): limit of -L_p/K = -(v_p(S_K)+gamma_p^out)/K from (59),(54),(58),(62),
    integrated exactly (rational arithmetic) over 1/3 < y < 2*lambda; comparison with Table 5 / (122).
 B. Inner function R(x) = -Gamma(x) - N(x) from (66)-(69), exact; comparison with (80)-(82), with rows of
    Tables 2-4, with (121) (integral over [3,20]); tail integral over [20,inf) (numerical) vs "A_0" (87).
 C. Direct evaluation, from definitions (49),(59),(60),(61),(62), of log m_{K,M}/K^2 at an admissible pair
    (K,M) = (320000,40) [K = 200 M^2], split by range of primes; comparison with the claimed A_M (89),
    with the corrected asymptotic value, and with the necessary condition for Prop. 5.1
        v_p(m_{K,M} S_K lead(Delta_K)) >= 0   (lead(Delta_K) from Lemma 2.3, eq. (10)).
 D. The leading-coefficient bound: log|lead F_K| / K^2 for various K, its limit, and C_* (95).
No claims are made by this file; it only computes.     Usage: python3 anand_zeta7_constants.py [A|B|C|D ...]
"""
import sys, math
from fractions import Fraction as Fr
from math import floor

AL, LAM, HH = Fr(3, 40), Fr(37, 40), Fr(49, 40)


def pos(x):
    return x if x > 0 else Fr(0)


# ------------------------------------------------------------------ Part A
def minus_gamma_out(y):
    """lim -gamma_p^out / K as p/K -> y, from (59) with (54),(58).  (gamma_out := 0 for p > K.)"""
    if y >= 1:
        return Fr(0)
    r = min(LAM, pos(1 + 6 * AL - 3 * y))
    if y > Fr(1, 2):
        v = 1 - y
        u = pos(AL + v - y)
        t = min(AL, v) + u
        return 9 * (1 - y) - 8 * t + min(r, y - AL + u)
    v = 1 - 2 * y
    u = pos(AL + v - y)
    t = min(AL, v) + u
    return 9 * (1 - y) - 16 * AL - 7 * t + min(r, y + u)


def S_part(y):
    """lim -v_p(S_K)/K, from (62):  -2 lam floor(1/y) + sum_j (2 lam - j y)_+ ."""
    return -2 * LAM * floor(1 / y) + sum(pos(2 * LAM - j * y) for j in range(1, 8))


def d_paper(y):
    """(72)"""
    if Fr(1, 3) < y < Fr(1, 2):
        return pos(1 + 6 * AL - 3 * y - pos(1 + 2 * AL - 3 * y))
    return Fr(0)


def T_def(y):
    return minus_gamma_out(y) + S_part(y)


def affine_pieces(f, lo, hi, cands):
    """exact piecewise-affine decomposition of f on [lo,hi] given candidate breakpoints; verified by a 3rd point."""
    pts = sorted(set([lo, hi] + [c for c in cands if lo < c < hi]))
    out = []
    for l, r in zip(pts[:-1], pts[1:]):
        x1, x2, x3 = (3 * l + r) / 4, (l + r) / 2, (l + 3 * r) / 4
        c = (f(x3) - f(x1)) / (x3 - x1)
        b = f(x1) - c * x1
        assert f(x2) == b + c * x2, ("not affine on", l, r)
        x4 = l + (r - l) / 97
        x5 = r - (r - l) / 89
        assert f(x4) == b + c * x4 and f(x5) == b + c * x5, ("not affine on", l, r)
        out.append((l, r, b, c))
    # merge
    merged = []
    for pc in out:
        if merged and merged[-1][2] == pc[2] and merged[-1][3] == pc[3]:
            merged[-1] = (merged[-1][0], pc[1], pc[2], pc[3])
        else:
            merged.append(pc)
    return merged


def partA():
    print("=== Part A: the outer integral I_out ===")
    # candidate breakpoints: all y = (a + b*alpha + c*lam)/k for small integers
    cands = set()
    for k in range(1, 9):
        for a in range(-3, 4):
            for b in range(-8, 9):
                for c in range(-2, 5):
                    y = (a + b * AL + c * LAM) / k
                    if Fr(1, 3) < y < 2 * LAM:
                        cands.add(y)
    lo, hi = Fr(1, 3), 2 * LAM
    pcs = affine_pieces(T_def, lo, hi, cands)
    I_def = sum(b * (r - l) + c * (r * r - l * l) / 2 for (l, r, b, c) in pcs)
    print("  T_def(y) := lim -(v_p(S_K)+gamma_p^out)/K  [(59),(62); this is the integrand that (60)-(61) actually produce]")
    print("  affine pieces of T_def (exact):   interval | b | c      (T = b + c y)")
    for (l, r, b, c) in pcs:
        print(f"    [{str(l):>7s}, {str(r):>7s}]  b = {str(b):>8s}  c = {str(c):>4s}")
    print(f"  I_out recomputed from (59),(62)       = {I_def} = {float(I_def):+.9f}")
    Id = sum(b * (r - l) + c * (r * r - l * l) / 2 for (l, r, b, c) in affine_pieces(d_paper, lo, hi, cands))
    print(f"  integral of d(y) of (72)               = {Id} = {float(Id):+.9f}")
    print(f"  variant 'R_0 := lim -gamma/K, then subtract d once more' (literal (71)+(73)): {float(I_def - Id):+.9f}")
    # Table 5 as printed
    T5 = [(Fr(1, 3), Fr(37, 100), Fr(-75, 8), 8), (Fr(37, 100), Fr(37, 80), Fr(-99, 8), 16),
          (Fr(37, 80), Fr(29, 60), Fr(-569, 40), 20), (Fr(29, 60), Fr(1, 2), Fr(-569, 40), 20),
          (Fr(1, 2), Fr(37, 60), Fr(-107, 10), 9), (Fr(37, 60), Fr(37, 40), Fr(-251, 20), 12),
          (Fr(37, 40), Fr(77, 80), Fr(-72, 5), 14), (Fr(77, 80), Fr(1), Fr(1), -2),
          (Fr(1), Fr(37, 20), Fr(-239, 20), -2)]
    I5 = sum(b * (r - l) + Fr(c) * (r * r - l * l) / 2 for (l, r, b, c) in T5)
    print(f"  Table 5 as printed, integrated          = {I5} = {float(I5):+.9f}")
    print(f"  preprint (122) claims I_out             = -11002997/720000 = {-11002997 / 720000:+.9f};  equal to Table-5 integral: {I5 == Fr(-11002997, 720000)}")
    print("  comparison on the Table-5 intervals (value at the midpoint):  interval | Table 5 | T_def | T_def - d")
    for (l, r, b, c) in T5:
        ym = (l + r) / 2
        print(f"    [{str(l):>7s}, {str(r):>7s}]  {float(b + c * ym):+9.4f}  {float(T_def(ym)):+9.4f}  {float(T_def(ym) - d_paper(ym)):+9.4f}")
    print("  on [1, 37/20]: gamma_out = 0 by definition (text before (60)), so T = (2 lam - y)_+ >= 0; Table 5 has -239/20 - 2y.")
    c1 = sum(b * (r - l) + Fr(c) * (r * r - l * l) / 2 for (l, r, b, c) in T5[-1:])
    print(f"    contribution of that single row to (122): {float(c1):+.6f};   correct contribution: {float((2 * LAM - 1) ** 2 / 2):+.6f}")
    return I_def


# ------------------------------------------------------------------ Part B
def ell(x, z):
    return floor(x - z) + floor(x + z) + 1


def dfun(u):
    return min(u, 1 - u)


def Gamma(x):
    """(66), exact for rational x."""
    T = floor(2 * HH * x)
    s = HH * x - Fr(T, 2)
    q = floor(2 * x)
    npl = (2 * x - q) / 2
    f, g = x - floor(x), AL * x - floor(AL * x)
    bps = sorted(set([Fr(0), Fr(1, 2)] + [z for z in (dfun(f), dfun(g)) if 0 < z < Fr(1, 2)]))
    integ = Fr(0)
    for l, r in zip(bps[:-1], bps[1:]):
        z = (l + r) / 2
        b = 4 * ell(AL * x, z)
        integ += (r - l) * (T - b) * (T + b - ell(x, z) - 7)
    return integ + s * (2 * T - q - 7) + pos(s - npl)


def Jfun(u):
    m = floor(2 * u)
    return m * u - Fr(m * (m + 1), 4)


def Nfun(x):
    return 2 * LAM * x * floor(x) - 16 * LAM * x * floor(AL * x) - 2 * Jfun(LAM * x)


def R(x):
    return -Gamma(x) - Nfun(x)


def R_via_80(x, cB=25):
    """(80)-(82); the printed coefficient of int B^2 is cB = 25"""
    fr = lambda u: u - floor(u)
    f, g, tau, sig, eta = fr(x), fr(AL * x), fr(2 * HH * x), fr(2 * x), fr(2 * LAM * x)
    eps = lambda u: 1 if u < Fr(1, 2) else -1
    df, dg = dfun(f), dfun(g)
    B2 = dg * (1 - 2 * dg)
    AB = eps(f) * eps(g) * (min(df, dg) - 2 * df * dg)
    F = 6 * LAM + 2 * LAM * f - 16 * LAM * g
    Q = (tau * (tau - sig) - pos(tau - sig) + eta * (1 - eta)) / 2 + cB * B2 - 4 * AB
    return x * F + Q


def R_breaks(lo, hi):
    cs = [Fr(1), Fr(2), AL, 2 * AL, 4 * AL, 8 * AL, LAM, 2 * LAM, HH, 2 * HH, 1 - AL, 1 + AL, 2 * (1 - AL),
          2 * (1 + AL), 2 * (1 + 2 * AL), 2 * (HH - 1), HH - 1, 2 * HH - 2 * LAM, 2 * HH - 1, 2 * HH + 1, 2 * HH - 2,
          2 * HH - AL, 2 * HH + AL, 1 - 2 * AL, 1 + 2 * AL]
    pts = set()
    for c in cs:
        k = floor(c * lo)
        while Fr(k) / c <= hi:
            if lo < Fr(k) / c < hi:
                pts.add(Fr(k) / c)
            k += 1
    return pts


def int_R_over_x3(lo, hi):
    pcs = affine_pieces(R, lo, hi, R_breaks(lo, hi))
    tot = Fr(0)
    for (l, r, b, c) in pcs:       # R = b + c x ;  int (b + c x)/x^3 = c(1/l - 1/r) + b/2 (1/l^2 - 1/r^2)
        tot += c * (1 / l - 1 / r) + b / 2 * (1 / l ** 2 - 1 / r ** 2)
    return tot, pcs


def partB():
    print("=== Part B: the inner function R(x) = -Gamma(x) - N(x) of (66)-(69) ===")
    bad = bad16 = 0
    import random
    random.seed(1)
    qmin, qmax = 10.0, -10.0
    for _ in range(4000):
        x = Fr(random.randint(3 * 10 ** 6, 60 * 10 ** 6), 10 ** 6)
        if R(x) != R_via_80(x):
            bad += 1
        if R(x) != R_via_80(x, 16):
            bad16 += 1
        f, g = x - floor(x), AL * x - floor(AL * x)
        qv = float(R(x) - x * (6 * LAM + 2 * LAM * f - 16 * LAM * g))
        qmin, qmax = min(qmin, qv), max(qmax, qv)
    print(f"  identity (80)-(82)  R(x) = x F(x) + Q(x) as printed (25 int B^2): mismatches at 4000 random rational x in [3,60]: {bad}")
    print(f"  same with 16 int B^2 in place of 25 int B^2 in (82):               mismatches: {bad16}")
    print(f"  observed range of Q = R - xF on the sample: [{qmin:+.4f}, {qmax:+.4f}]   [(83) claims -1/2 <= Q <= 27/8: still true]")
    I320, pcs = int_R_over_x3(Fr(3), Fr(20))
    print(f"  number of affine pieces of R on [3,20] (after merging): {len(pcs)}   [preprint: 192 intervals]")
    print(f"  int_3^20 R(x)/x^3 dx (exact)  = {I320}\n                               = {float(I320):+.12f}")
    claim = Fr(120667538150827615739401, 708353029541942104320000)
    print(f"  preprint (121)                = {float(claim):+.12f};   equal: {I320 == claim}")
    # table rows
    rows = [(Fr(3), Fr(70, 23), Fr(51, 20), 1), (Fr(120, 37), Fr(140, 43), Fr(81, 10), -17), (Fr(160, 49), Fr(10, 3), Fr(55, 8), -13),
            (Fr(4), Fr(200, 49), Fr(67, 20), -7), (Fr(5), Fr(220, 43), Fr(27, 10), -12), (Fr(12), Fr(520, 43), Fr(9, 2), None)]
    print("  spot check of printed table rows (a_i, b_i) against R:   interval | printed a_i x + b_i at midpoint | R(midpoint) | slope of R")
    for (l, r, a, b) in rows:
        xm = (l + r) / 2
        x1, x2 = (2 * l + r) / 3, (l + 2 * r) / 3
        slope = (R(x2) - R(x1)) / (x2 - x1)
        pr = '   n/a  ' if b is None else f"{float(a * xm + b):+8.4f}"
        print(f"    [{str(l):>7s},{str(r):>7s}]  printed a={str(a):>6s} b={str(b):>4s}: {pr}   R = {float(R(xm)):+8.4f}   slope = {str(slope)}")
    # tail
    X = Fr(20)
    tot = 0.0
    step = 40
    Xmax = 20 + 40 * 60
    while X < Xmax:
        v, _ = int_R_over_x3(X, X + step)
        tot += float(v)
        X += step
    # beyond Xmax: R = xF + Q, mean F = -lam; error O(1/X^2)
    tail = tot - float(LAM) / Xmax
    print(f"  int_20^{Xmax} R/x^3 dx = {tot:+.8f};  remaining tail ~ -lam/X = {-float(LAM) / Xmax:+.8f} (error O(30/X^2) ~ {30 / Xmax ** 2:.1e})")
    print(f"  int_20^inf R(x)/x^3 dx ~ {tail:+.6f}")
    A0 = Fr(9928298118277006344769, 7535670527041937280000)
    print(f"  preprint: 'A_0 = int_20^inf R/x^3 dx' with A_0 (87) = {float(A0):+.10f}")
    print("            [(87) is digit-for-digit Fauzan's zeta(5) constant A_* (his (5.19)) = his I_out + his int_3^20 - 2689/48000]")
    chk = Fr(127751, 96000) + Fr(322437603634266857629, 7535670527041937280000) - Fr(2689, 48000)
    print(f"            check: 127751/96000 + 322437603634266857629/7535670527041937280000 - 2689/48000 == A_0 : {chk == A0}")
    return float(I320), tail


# ------------------------------------------------------------------ Part C
def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(n + 1) if sieve[i]]


def legendre(m, p):
    v = 0
    while m:
        m //= p
        v += m
    return v


def sum_floor_2i(h, q):
    """sum_{i=1}^{h-1} floor(2i/q)"""
    tot = 0
    j = 1
    while True:
        c = h - (j * q + 1) // 2      # #{i<=h-1 : 2i >= j q} = h-1 - ceil(jq/2) + 1
        if c <= 0:
            break
        tot += c
        j += 1
    return tot


def vp_SK(K, N, h, p):
    v = 2 * h * legendre(K, p) - 16 * h * legendre(N, p)
    q = p
    while q <= 2 * h:
        v -= 2 * sum_floor_2i(h, q)
        q *= p
    if p == 2:
        v += 2 * (h - 1)
    return v


def Sf(X, p):
    """sum_{n=0}^{X} v_p(n!)"""
    if X < 0:
        return 0
    tot = 0
    q = p
    while q <= X:
        Qf = X // q
        tot += q * Qf * (Qf - 1) // 2 + Qf * (X - q * Qf + 1)
        q *= p
    return tot


def vp_lead_Delta(K, N, p):
    """v_p of prod_{j=N+1}^K j^6 D_N(-j^2)^7,  D_N(-j^2) = +-[(j-1)!/(j-N-1)!][(j+N)!/j!]."""
    return 6 * (legendre(K, p) - legendre(N, p)) + 7 * ((Sf(K - 1, p) - Sf(N - 1, p)) - Sf(K - N - 1, p)
                                                         + (Sf(K + N, p) - Sf(2 * N, p)) - (Sf(K, p) - Sf(N, p)))


def gamma_out(K, N, h, p):
    if not (p >= 11 and 3 * p > K and p <= K and p * p > 2 * K and 7 * N <= 2 * p - 2):
        return None
    v = K - p * (K // p)
    u = max(0, N + v - p + 1)
    tp = min(N, v) + u
    rp = min(h, max(0, K + 6 * N - 3 * p + 4))
    if K < 2 * p:
        return -9 * (K - p) + 8 * tp - 1 - min(rp, p - 1 - N + u)
    return -9 * (K - p) + 3 + 16 * N + 7 * tp - min(rp, p + u)


def gamma_in(K, N, h, M, p):
    """(49) with Definitions 4.1-4.3; returns (gamma, min L_a)."""
    import numpy as np
    m = (p - 1) // 2
    mN, mK = N // p, K // p
    L0 = 6 * M + 10
    a = np.arange(1, m + 1, dtype=np.int64)
    lK = (K - a) // p + 1 + (K + a) // p
    lN = (N - a) // p + 1 + (N + a) // p
    assert int(lN.sum()) == N - mN and int(lK.sum()) == K - mK
    tot = h - L0 + 4 * (N - mN)
    T, E = divmod(tot, m)
    order = np.argsort(-lK, kind='stable')
    eps = np.zeros(m, dtype=np.int64)
    eps[order[:E]] = 1
    La = T - 4 * lN + eps
    Za = T + eps
    assert int(La.sum()) + L0 == h
    # 2*sum_i w_{a,i} = La(La-1) + La(8 lN - lK - 6)
    two_sum = int((La * (La - 1) + La * (8 * lN - lK - 6)).sum())
    # zero class: w_{0,i} = min(2i + 8 mN - mK + 1/2, min_c [Z_c - (lK(c)+6)/2]);  work with doubled values
    cap2 = int((2 * Za - lK - 6).min())
    z2 = 0
    for i in range(L0):
        z2 += min(4 * i + 16 * mN - 2 * mK + 1, cap2)
    return two_sum + z2, int(La.min())


def ilog(p, x):
    k = 0
    while p ** (k + 1) <= x:
        k += 1
    return k


def partC(K=320000, M=40, I_out=None, I320=None, tail=None):
    N, h = 3 * K // 40, 37 * K // 40
    assert K % 40 == 0 and M % 40 == 0 and K >= 200 * M * M
    print(f"=== Part C: m_(K,M) from the definitions (60),(61) at K={K}, M={M} (admissible: K >= 200 M^2 = {200 * M * M}) ===")
    P = primes_upto(2 * h)
    parts = {'crude  pM<=K': 0.0, 'inner  K/M<p<=K/3': 0.0, 'outer  K/3<p<=K': 0.0, 'outer  K<p<=2h': 0.0}
    viol = []
    minLa = 10 ** 9
    chk70 = []
    for p in P:
        lp = math.log(p)
        vS = vp_SK(K, N, h, p)
        if p * M <= K:
            L = -8 * h * ilog(p, 5 * K) - h * {2: 4, 3: 2, 5: 1}.get(p, 0)
            key = 'crude  pM<=K'
        elif 3 * p <= K:
            g, mla = gamma_in(K, N, h, M, p)
            minLa = min(minLa, mla)
            L = vS + g
            key = 'inner  K/M<p<=K/3'
            if len(chk70) < 6 and p > len(chk70) * 18000 + 8000:
                x = Fr(K, p)
                chk70.append((p, g, float(p * Gamma(x)), vS, float(p * Nfun(x))))
        elif p <= K:
            g = gamma_out(K, N, h, p)
            assert g is not None
            L = vS + g
            key = 'outer  K/3<p<=K'
        else:
            L = vS
            key = 'outer  K<p<=2h'
        parts[key] += -L * lp
        vl = vS + vp_lead_Delta(K, N, p)
        if L > vl:
            viol.append((p, key.split()[0], L, vl))
    K2 = K * K
    tot = sum(parts.values())
    print("  contribution to log m_(K,M) / K^2  ( = - sum_p L_p log p / K^2 ):")
    for k, v in parts.items():
        print(f"    {k:20s}: {v / K2:+.6f}")
    print(f"    {'TOTAL':20s}: {tot / K2:+.6f}        <-- log m_(K,M)/K^2 computed from the preprint's own definitions")
    lamf = float(LAM)
    B0 = float(Fr(9928298118277006344769, 7535670527041937280000)) - 11002997 / 720000 + float(Fr(120667538150827615739401, 708353029541942104320000))
    AM = B0 + 9 * lamf / M - 11549 / (720 * M * M) + 44 / M ** 3
    print(f"  preprint's claimed bound A_M (89) at M={M}: {AM:+.6f}    (A_200 = -13.752854, A_100000 = -13.794001)")
    print(f"  min over inner primes and classes of L_a (must be >= 1): {minLa}")
    print("  check of (70) at a few inner primes:  p | gamma_in (49) | p*Gamma(K/p) | v_p(S_K) | p*N(K/p)")
    for (p, g, pg, vS, pn) in chk70:
        print(f"    {p:7d} | {g:10d} | {pg:13.2f} | {vS:10d} | {pn:13.2f}")
    if I_out is not None:
        # asymptotic prediction from the proof of Prop 5.4 with recomputed integrals
        Iin, _ = int_R_over_x3(Fr(3), Fr(M)) if M <= 60 else (None, None)
        if Iin is not None:
            pred = 8 * lamf / M + float(Iin) + float(I_out)
            print(f"  asymptotic value by the preprint's own argument with recomputed integrals: 8 lam/M + int_3^M R/x^3 + I_out = "
                  f"{8 * lamf / M:+.6f} {float(Iin):+.6f} {float(I_out):+.6f} = {pred:+.6f}")
        I200, _ = int_R_over_x3(Fr(3), Fr(200))
        print(f"  corrected constant for M=200 (same formula): 8 lam/200 + int_3^200 R/x^3 + I_out = {8 * lamf / 200:+.6f} {float(I200):+.6f} {float(I_out):+.6f}"
              f" = {8 * lamf / 200 + float(I200) + float(I_out):+.6f}     [preprint: A_200 = -13.752854]")
        if I320 is not None and tail is not None:
            print(f"  corrected constant for M -> infinity: int_3^inf R/x^3 + I_out = {I320 + tail:+.6f} {float(I_out):+.6f} = {I320 + tail + float(I_out):+.6f}"
                  f"     [preprint: B_0 = -13.794084]")
            print(f"  with the preprint's real constant U = +1.76: corrected A_200 + U = {8 * lamf / 200 + float(I200) + float(I_out) + 1.76:+.4f} > 0;"
                  f"  with the numerical equilibrium value -0.7619: {8 * lamf / 200 + float(I200) + float(I_out) - 0.7619:+.4f} > 0")
    print(f"  necessary condition for Prop. 5.1 (Q in Z[X]):  L_p <= v_p(S_K) + v_p(lead Delta_K) for every p <= 2h")
    print(f"    number of primes violating it: {len(viol)}  (of {len(P)} primes <= 2h)")
    by = {}
    for (p, key, L, vl) in viol:
        by.setdefault(key, []).append((p, L, vl))
    for key, lst in by.items():
        worst = max(lst, key=lambda t: t[1] - t[2])
        print(f"    range {key}: {len(lst)} primes; first p={lst[0][0]} (L_p={lst[0][1]}, v_p(lead F_K)={lst[0][2]}); worst excess {worst[1] - worst[2]} at p={worst[0]}")
    ll = log_lead_F(K)
    print(f"  log|lead F_K|/K^2 = {ll / K2:+.6f}  => integrality of m*lead(F_K) forces log m/K^2 >= {-ll / K2:+.6f}")


# ------------------------------------------------------------------ Part D
def log_SK(K):
    N, h = 3 * K // 40, 37 * K // 40
    return (2 * h * math.lgamma(K + 1) + (h - 1) * math.log(4) - 16 * h * math.lgamma(N + 1)
            - 2 * sum(math.lgamma(2 * i + 1) for i in range(1, h)))


def log_lead_Delta(K):
    N, h = 3 * K // 40, 37 * K // 40
    lg = math.lgamma
    # sum_{j=N+1}^K [ 6 log j + 7 ( log (j-1)! - log (j-N-1)! + log (j+N)! - log j! ) ]
    tot = 6 * (lg(K + 1) - lg(N + 1))
    for j in range(N + 1, K + 1):
        tot += 7 * (lg(j) - lg(j - N) + lg(j + N + 1) - lg(j + 1))
    return tot


def log_lead_F(K):
    return log_SK(K) + log_lead_Delta(K)


def partD():
    import mpmath as mp
    print("=== Part D: the leading coefficient of F_K = S_K Delta_K ===")
    al, lam = mp.mpf(3) / 40, mp.mpf(37) / 40
    Cstar = -2 * lam + 16 * al * lam * (1 - mp.log(al)) + 3 * lam ** 2 - 2 * lam ** 2 * mp.log(2 * lam)
    dbl = mp.quad(lambda u: mp.quad(lambda v: mp.log(u * u - v * v), [0, al]), [al, 1])
    clead = Cstar + 7 * dbl
    print(f"  C_* (95) = {mp.nstr(Cstar, 15)}   [preprint: 3.649335128751705]")
    print(f"  lim log|lead F_K|/K^2 = C_* + 7 int_al^1 int_0^al log(u^2-v^2) dv du = {mp.nstr(clead, 10)}")
    print("     K   | log S_K/K^2 - (-14 al lam) log K | [paper (107) upper bound: C_* + (6h log K + 6h)/K^2] | log|lead F_K|/K^2")
    for K in (40, 80, 120, 160, 400, 4000, 40000, 320000):
        N, h = 3 * K // 40, 37 * K // 40
        lS = log_SK(K)
        a = lS / K ** 2 + 14 * float(al * lam) * math.log(K)
        b = float(Cstar) + (6 * h * math.log(K) + 6 * h) / K ** 2
        print(f"  {K:6d} | {a:+.6f} | {b:+.6f} | {log_lead_F(K) / K ** 2:+.6f}")
    print("  Consequence: if Q = m F_K has integer coefficients then |m lead(F_K)| >= 1, i.e. log m >= -log|lead F_K| ~ -2.88 K^2,")
    print("  whereas Prop. 5.4 claims limsup log m_(K,M)/K^2 <= A_M ~ -13.75.  Prop. 5.1 and Prop. 5.4 (with (89)) cannot both hold.")


if __name__ == '__main__':
    what = sys.argv[1:] or ['A', 'B', 'C', 'D']
    I_out = I320 = tail = None
    if 'A' in what or 'C' in what:
        I_out = partA()
    if 'B' in what:
        I320, tail = partB()
    if 'C' in what:
        partC(I_out=I_out, I320=I320, tail=tail)
    if 'D' in what:
        partD()
