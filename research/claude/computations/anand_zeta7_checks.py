#!/usr/bin/env python3
"""
anand_zeta7_checks.py -- numerical checks for the referee report on P. Anand, "zeta(7) is irrational".

 1. Prop. 2.4 (positive integral representation), 35 digits (mpmath):
      int_0^inf y^{2e} w_7(y) dy                 vs (5)
      int_0^inf w_7(y)/(y^2+j^2) dy              vs (6) [constant -1/6]  and vs the variant with -1/2
      int_0^inf D_N(y^2)^8 y^{2e}/D_K(y^2) w_7   vs alpha_e + beta_e zeta(7) from hankel_lib.moments_F1(7,K,N,8,h)
    (the last line ties the exact matrices used in anand_zeta7_exact.py to the positive weight).
 2. Real side: the comparison measure of Table 1 / Lemma 6.1: mass, I(rho) (110), C_* (95),
    sup_t [2U^rho - V] on a fine grid, the resulting constant; and the true constant from a numerical
    solution of the weighted equilibrium problem (hankel_hadamard.solve_equilibrium, a=8).
 3. Elementary lower bound (all steps are inequalities, see the report):
      Delta_K(zeta7)/|lead Delta_K| >= exp( -c(a,b) K^2 - h (log(2h) - 2 log 2 - 3 log a - log((2pi)^6/720)) ),
      c(a,b) = lam[ -8 al log a + log(1+b) + 2 pi sqrt b ] - lam^2 log((b-a)/4),
    compared with the claimed upper bound exp(-19000 n^2) = exp(-11.875 K^2) for EVERY n >= 1 and, where the
    exact value is available (n <= 4, from anand_zeta7_exact.py output), with the exact value.
No claims are made by this file; it only computes.
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
import mpmath as mp


def part1():
    mp.mp.dps = 40
    print("=== 1. Prop. 2.4: the functional at X = zeta(7) against the weight w_7 ===")
    twopi = 2 * mp.pi

    def w7(y):
        q = mp.exp(-twopi * y)
        return twopi ** 6 / 360 * y ** 7 * q * (1 + 57 * q + 302 * q ** 2 + 302 * q ** 3 + 57 * q ** 4 + q ** 5) / (1 - q) ** 7

    pts = [0, mp.mpf(1) / 8, mp.mpf(1) / 2, 1, 2, 4, 8, 16, 32, 64]
    worst = 0
    for e in range(0, 7):
        val = mp.quad(lambda y: y ** (2 * e) * w7(y), pts)
        ref = (-1) ** e * mp.bernoulli(2 * e + 2) * (2 * e + 3) * (2 * e + 4) * (2 * e + 5) * (2 * e + 6) * (2 * e + 7) / 720
        worst = max(worst, abs(val - ref) / abs(ref))
    print(f"  moments e=0..6: max relative error vs (5): {mp.nstr(worst, 3)}")
    z7 = mp.zeta(7)
    worst6 = 0
    for j in range(1, 9):
        val = mp.quad(lambda y: w7(y) / (y * y + j * j), pts)
        Hj = sum(mp.mpf(1) / v ** 7 for v in range(1, j + 1))
        ref6 = j ** 6 * (z7 - Hj) + mp.mpf(1) / (2 * j) - mp.mpf(1) / 6
        ref2 = j ** 6 * (z7 - Hj) + mp.mpf(1) / (2 * j) - mp.mpf(1) / 2
        worst6 = max(worst6, abs(val - ref6))
        if j <= 2:
            print(f"  j={j}: integral = {mp.nstr(val, 32)}\n        (6) with -1/6 = {mp.nstr(ref6, 32)}   diff {mp.nstr(val - ref6, 3)};   variant with -1/2: diff {mp.nstr(val - ref2, 6)}")
    print(f"  poles j=1..8: max |integral - (6)| = {mp.nstr(worst6, 3)}")
    from hankel_lib import moments_F1
    K, N, h = 5, 1, 4
    al, be = moments_F1(7, K, N, 8, h)
    worst = 0
    for e in range(2 * h - 1):
        def f(y, e=e):
            t = y * y
            num = mp.mpf(1)
            for j in range(1, N + 1):
                num *= (t + j * j)
            den = mp.mpf(1)
            for j in range(1, K + 1):
                den *= (t + j * j)
            return num ** 8 * t ** e / den * w7(y)
        val = mp.quad(f, pts)
        ref = mp.mpf(int(al[e].p)) / int(al[e].q) + z7 * mp.mpf(int(be[e].p)) / int(be[e].q)
        worst = max(worst, abs(val - ref) / abs(val))
    print(f"  Gram moments of D_N^8 t^e/D_K (K=5,N=1, e=0..6) vs exact alpha_e + beta_e zeta(7): max rel. error {mp.nstr(worst, 3)}")
    print("     (cancellation in alpha_e + beta_e zeta(7) limits this comparison; 40-digit arithmetic)")


def part2():
    import numpy as np
    print("=== 2. Real side: Lemma 6.1 / Table 1, and the true equilibrium constant ===")
    L = [0.87, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.25, 0.20, 0.15, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02]
    c = [0.10, 0.08, 0.06] + [0.05] * 11 + [0.015, 0.12]
    al, lam = 3 / 40, 37 / 40
    print(f"  mass of rho = {sum(c):.6f} (lam = {lam})")
    S = np.cumsum(c)
    I = sum((S[j] ** 2 - (S[j - 1] ** 2 if j else 0)) * math.log(L[j] / 4) for j in range(16))
    Cs = -2 * lam + 16 * al * lam * (1 - math.log(al)) + 3 * lam ** 2 - 2 * lam ** 2 * math.log(2 * lam)
    print(f"  I(rho) by (110) = {I:.12f}   [preprint: -3.143325089379538]")
    print(f"  C_* (95)        = {Cs:.12f}   [preprint:  3.649335128751705]")

    def U(t):
        tot = np.zeros_like(t)
        for Lj, cj in zip(L, c):
            inside = np.log(Lj / 4) * np.ones_like(t)
            with np.errstate(invalid='ignore', divide='ignore'):
                outside = np.log((t - Lj / 2 + np.sqrt(np.maximum(t * (t - Lj), 0))) / 2)
            tot += cj * np.where(t <= Lj, inside, outside)
        return tot

    def V(t):
        st = np.sqrt(t)
        J = lambda a: a * np.log(t + a * a) - 2 * a + 2 * st * np.arctan(a / st)
        return 2 * np.pi * st + J(1.0) - 8 * J(al)

    t = np.concatenate([np.linspace(1e-9, 0.05, 400001), np.linspace(0.05, 3, 1200001), np.linspace(3, 400, 400001)])
    g = 2 * U(t) - V(t)
    k = int(np.argmax(g))
    print(f"  sup_t [2U^rho - V] on a grid of {len(t)} points in (0,400]: {g[k]:.6f} at t = {t[k]:.5f}   [preprint: about -5.4719 at t ~ 0.2648; M_0 = -5.45]")
    print(f"  lam*M_0 - I(rho) + C_* with M_0 = -5.45: {lam * (-5.45) - I + Cs:.9f}   [preprint: 1.751410218131243 <= U = 1.76]")
    print(f"  same with the grid supremum in place of M_0: {lam * g[k] - I + Cs:.6f}")
    from hankel_hadamard import solve_equilibrium
    for n in (800, 1600):
        E, m, edges, F = solve_equilibrium(lam, a=8, alpha=al, Y=2.5, n=n)
        sup = edges[1:][m > 1e-12].max()
        print(f"  equilibrium problem (a=8, al=3/40, lam=37/40), {n} cells: E* = {E:.6f}, Robin const = {F:.5f}, support t <= {sup ** 2:.4f};"
              f"  true real constant C_* - E* = {Cs - E:+.6f}")
    print("  exact values log F_K(zeta7)/K^2 (anand_zeta7_exact.py): K=40: -0.655371, 80: -0.707302, 120: -0.725215, 160: -0.734284")


def part3():
    print("=== 3. Elementary lower bound for Delta_K(zeta7)/|lead Delta_K| (hence for P_K(zeta7)) ===")
    al, lam = 3 / 40, 37 / 40

    def cG(a, b):
        return lam * (-8 * al * math.log(a) + math.log(1 + b) + 2 * math.pi * math.sqrt(b)) - lam ** 2 * math.log((b - a) / 4)
    for (a, b) in ((1 / 18, 1 / 3), (3 / 40, 1 / 5)):
        print(f"  c(a,b) at (a,b)=({a:.5f},{b:.5f}): {cG(a, b):.6f}")
    a, b = 3 / 40, 1 / 5
    c = cG(a, b)
    exact = {1: (-0.655371, 2.944310, 1.451854), 2: (-0.707302, 2.907567, 1.578523), 3: (-0.725215, 2.894645, 1.609096), 4: (-0.734284, 2.888061, 1.619459)}
    print("   n |    K | lower bound for log P_K(zeta7) | claimed upper bound -19000 n^2 | exact log(Delta(zeta7)/|lead|) | exact log P_K(zeta7)")
    ok = True
    for n in list(range(1, 11)) + [100, 1000, 10 ** 4, 10 ** 6]:
        K, h = 40 * n, 37 * n
        lb = -c * K * K + h * (math.log((2 * math.pi) ** 6 / 720) + 3 * math.log(a) + 2 * math.log(2) - math.log(2 * h))
        cl = -19000 * n * n
        ok = ok and lb > cl
        ex = exact.get(n)
        e1 = f"{(ex[0] - ex[1]) * K * K:14.1f}" if ex else "       .      "
        e2 = f"{(ex[0] + ex[2]) * K * K:+12.1f}" if ex else "      .     "
        if n <= 6 or n >= 100:
            print(f"  {n:7d} | {K:8d} | {lb:16.1f} | {cl:16.1f} | {e1} | {e2}")
    # for all n: lb - cl = (11.875 - c) K^2 - h(...) ; positivity for all n >= 1:
    worst = min(((-c * (40 * n) ** 2 + 37 * n * (math.log((2 * math.pi) ** 6 / 720) + 3 * math.log(a) + 2 * math.log(2) - math.log(74 * n))) + 19000 * n * n) / n ** 2 for n in range(1, 200000))
    print(f"  min over 1 <= n < 200000 of (lower bound + 19000 n^2)/n^2 = {worst:.2f} > 0;  for larger n the difference is >= (11.875 - {c:.4f}) K^2 - K(2 + log(2K)) > 0")
    print(f"  => for every n >= 1 the lower bound exceeds the claimed upper bound: {ok}")


def part3b(n=1):
    """step-by-step check of the bound at K=40: d_k = det G_{k+1}/det G_k >= rho_min (B-A)^{2k+1} (k!)^4/((2k)!^2 (2k+1))."""
    from hankel_lib import moments_F1
    from flint import arb, arb_mat, ctx
    K, N, h = 40 * n, 3 * n, 37 * n
    al, be = moments_F1(7, K, N, 8, h)
    ctx.prec = 4000 + 80 * K * K
    z = arb(7).zeta()
    mom = [arb(a.p) / arb(a.q) + z * arb(b.p) / arb(b.q) for a, b in zip(al, be)]
    a_, b_ = 3 / 40, 1 / 5
    A, B = a_ * K * K, b_ * K * K
    log_rho = (8 * N * math.log(A) - 2 * K * math.log(K) - K * math.log(1 + b_) + math.log((2 * math.pi) ** 6 / 720)
               + 3 * math.log(A) - 2 * math.pi * K * math.sqrt(b_))
    prev = arb(1)
    worst, ok = 1e99, True
    for k in range(h):
        M = arb_mat(k + 1, k + 1, [mom[i + j] for i in range(k + 1) for j in range(k + 1)])
        d = M.det()
        ok = ok and bool(d > 0)
        ldk = float((d / prev).log().mid())
        prev = d
        lb = log_rho + (2 * k + 1) * math.log(B - A) + 4 * math.lgamma(k + 1) - 2 * math.lgamma(2 * k + 1) - math.log(2 * k + 1)
        ok = ok and ldk >= lb
        worst = min(worst, ldk - lb)
    print(f"  3b. K={K}: all leading principal minors of G_K(zeta7) positive and every d_k >= its lower bound: {ok};"
          f" min slack in log d_k: {worst:.2f};  log det G_K(zeta7) = {float(prev.log().mid()):.6f}")


if __name__ == '__main__':
    what = sys.argv[1:] or ['1', '2', '3']
    if '1' in what:
        part1()
    if '2' in what:
        part2()
    if '3' in what:
        part3()
        part3b(1)
