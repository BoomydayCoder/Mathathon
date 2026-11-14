#!/usr/bin/env python3
"""
hankel_hadamard.py -- analytic side for s = 3, W = 1, h = K.

Delta_K(zeta(3)) = det of the Gram matrix of 1, t, ..., t^{K-1} for the positive weight
    omega_K(y) dy = w_3(y)/D_K(y^2) dy on (0,inf),  t = y^2,  w_3(y) = pi^2 y^3 cosh(pi y)/sinh(pi y)^3.
Hadamard's inequality for Gram matrices: for ANY monic p_n of degree n,
    Delta_K(zeta(3)) <= prod_{n<K} ||p_n||^2,   ||p||^2 = int p(y^2)^2 omega_K(y) dy = mu_{zeta(3)}(p^2/D_K).
We take p_n(t) = prod_{j=1}^{n} (t - (theta j)^2).

Part 1: the one-variable inequality (H_theta)
    2 phi(r) <= 2 pi theta r + log(1+theta^2 r^2) + 2 theta r arctan(1/(theta r)),  r > 0,
    phi(r) = (1+r)log(1+r) + (1-r) log|1-r|,
  which (in the scaling limit) says that p_n(y^2)^2 omega_K(y) is maximal at y = O(1).  Margin on a grid.
Part 2: exact finite-K check: sum_n log ||p_n||^2 (computed exactly as a_n + b_n zeta(3), then evaluated
  in high precision) against log Delta_K(zeta(3)) and against the predicted (2 log theta - 1) K^2.
Part 3: the equilibrium problem: closed-form endpoint b(lambda): 3b + 1 - sqrt(1+b^2) = 2 lambda,
  and numerical energies E*(lambda) (needs energy solver below); E*(1) vs (9/2) log 2.
Usage: python3 hankel_hadamard.py [Kmax]
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
import numpy as np
from hankel_lib import *
from fractions import Fraction as Fr


def H_margin(theta, rmax=60.0, n=2000001):
    r = np.linspace(1e-9, rmax, n)
    with np.errstate(divide='ignore', invalid='ignore'):
        phi = (1 + r) * np.log(1 + r) + np.where(np.abs(1 - r) > 0, (1 - r) * np.log(np.abs(1 - r)), 0.0)
    rhs = 2 * np.pi * theta * r + np.log(1 + theta ** 2 * r ** 2) + 2 * theta * r * np.arctan2(1.0, theta * r)
    d = rhs - 2 * phi
    # margin relative to r (we need rhs - 2 phi >= delta * theta * r for some delta > 0, r >= r0)
    k = np.argmin(d[1000:]) + 1000
    rel = d[1000:] / (theta * r[1000:])
    k2 = np.argmin(rel) + 1000
    return d[k], r[k], rel.min(), r[k2]


def part1():
    print("Part 1: inequality (H_theta); c(theta) = 1 - 2 log(theta)")
    print("  theta | c(theta) | min_r [RHS - LHS] at r | min_r [RHS-LHS]/(theta r) (r>=0.03) at r")
    for theta in (0.38, 0.39, 0.40, 0.42, 0.45, 0.5):
        m, rm, rel, rr = H_margin(theta)
        print(f"  {theta:.2f}  | {1 - 2 * math.log(theta):.4f}  | {m:+.5f} at r={rm:.3f} | {rel:+.5f} at r={rr:.3f}")


def norms_exact(K, theta_num, theta_den, s=3):
    """returns list of (a_n, b_n) with ||p_n||^2 = a_n + b_n zeta(3), p_n = prod_{j<=n}(t - (theta j)^2)."""
    al, be = moments_F1(s, K, 0, 0, K)
    out = []
    P = fmpq_poly([1])
    th2 = fmpq(theta_num * theta_num, theta_den * theta_den)
    for n in range(K):
        if n >= 1:
            P = P * fmpq_poly([-th2 * n * n, 1])
        sq = (P * P).coeffs()
        a = sum((sq[e] * al[e] for e in range(len(sq))), fmpq(0))
        b = sum((sq[e] * be[e] for e in range(len(sq))), fmpq(0))
        out.append((a, b))
    return out, al, be


def part2(Ks, thetas):
    print("Part 2: exact finite-K Hadamard check (s=3, W=1, h=K)")
    print("  K | log Delta_K(zeta3)/K^2 | theta | sum_n log||p_n||^2 / K^2 | predicted limit 2 log(theta) - 1 | Hadamard holds")
    for K in Ks:
        al, be = moments_F1(3, K, 0, 0, K)
        cs = det_poly_affine(hankel(al, K), hankel(be, K))
        prec = 6000 + 80 * K * K
        val = eval_poly_at_zeta(cs, 3, prec)
        ld, _ = log_abs_arb(val)
        for (tn, td) in thetas:
            nr, _, _ = norms_exact(K, tn, td)
            old = ctx.prec
            ctx.prec = prec
            z = arb(3).zeta()
            tot = 0.0
            ok = True
            for (a, b) in nr:
                v = arb(a.p) / arb(a.q) + z * arb(b.p) / arb(b.q)
                if not (v > 0):
                    ok = False
                tot += float(v.log().mid())
            ctx.prec = old
            print(f"  {K:3d} | {ld / K ** 2:9.4f} | {tn}/{td} | {tot / K ** 2:9.4f} | {2 * math.log(tn / td) - 1:8.4f} | {ok and tot >= ld}")
        sys.stdout.flush()


# ---------- equilibrium problem (numerical) ----------
def Fint(x):
    ax = np.abs(x)
    with np.errstate(divide='ignore', invalid='ignore'):
        r = 0.5 * x * x * np.log(ax) - 0.75 * x * x
    return np.where(ax == 0, 0.0, r)


def cell_log_matrix(edges):
    a = edges[:-1]; b = edges[1:]
    A1, A2 = np.meshgrid(a, a, indexing='ij'); B1, B2 = np.meshgrid(b, b, indexing='ij')
    Mminus = -(Fint(B1 - B2) - Fint(B1 - A2) - Fint(A1 - B2) + Fint(A1 - A2))
    Mplus = (Fint(B1 + B2) - Fint(B1 + A2) - Fint(A1 + B2) + Fint(A1 + A2))
    return Mminus + Mplus


def Vfun(y, a=0.0, alpha=0.0):
    t = y * y
    def J(c):
        st = np.sqrt(t)
        with np.errstate(divide='ignore', invalid='ignore'):
            r = c * np.log(t + c * c) - 2 * c + 2 * st * np.arctan(c / np.where(st == 0, 1, st))
        return np.where(st == 0, c * np.log(c * c) - 2 * c, r)
    v = 2 * np.pi * y + J(1.0)
    if a != 0:
        v = v - a * J(alpha)
    return v


def solve_equilibrium(lam, a=0.0, alpha=0.0, Y=2.5, n=800):
    """min over rho>=0 of mass lam on y in [0,Y] of int V drho - int int log|y^2-z^2| drho drho.
       returns (E*, masses, edges, Robin constant F = V - 2U on the support)."""
    s_ = np.linspace(0, 1, n + 1)
    edges = Y * s_ ** 1.5
    mid = 0.5 * (edges[:-1] + edges[1:]); wid = np.diff(edges)
    M = cell_log_matrix(edges)
    gx, gw = np.polynomial.legendre.leggauss(6)
    Vc = np.zeros(n)
    for k in range(6):
        yy = mid + 0.5 * wid * gx[k]
        Vc += 0.5 * wid * gw[k] * Vfun(yy, a, alpha)
    Q = -M / np.outer(wid, wid)
    c = Vc / wid
    active = np.ones(n, bool)
    m = np.full(n, lam / n); mu = 0.0
    for it in range(400):
        idx = np.where(active)[0]
        k = len(idx)
        KK = np.zeros((k + 1, k + 1)); rhs = np.zeros(k + 1)
        KK[:k, :k] = 2 * Q[np.ix_(idx, idx)]; KK[:k, k] = 1; KK[k, :k] = 1
        rhs[:k] = -c[idx]; rhs[k] = lam
        sol = np.linalg.solve(KK, rhs)
        mm = np.zeros(n); mm[idx] = sol[:k]; mu = sol[k]
        neg = mm < -1e-13
        if neg.any():
            active = active & ~neg
            continue
        grad = 2 * Q @ mm + c + mu
        viol = (~active) & (grad < -1e-10)
        if viol.any():
            active = active | viol
            continue
        m = mm
        break
    E = m @ Q @ m + c @ m
    return E, m, edges, -mu


def part3():
    print("Part 3: equilibrium problem  E*(lam) = min [ int V drho - I(rho) ],  V(eta) = 2 pi eta + int_0^1 log(eta^2+u^2) du")
    print("  lam | E*(lam) numeric | right end of support (numeric) | closed form b: 3b+1-sqrt(1+b^2) = 2 lam")
    for lam in (0.5, 0.75, 1.0, 1.25):
        E, m, edges, F = solve_equilibrium(lam, n=1000)
        bnum = edges[1:][m > 1e-12].max()
        # solve 3b + 1 - sqrt(1+b^2) = 2 lam
        lo, hi = 0.0, 10.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if 3 * mid + 1 - math.sqrt(1 + mid * mid) < 2 * lam:
                lo = mid
            else:
                hi = mid
        print(f"  {lam:4.2f} | {E:.6f} | {bnum:.4f} | {lo:.4f}")
    print(f"  (9/2) log 2 = {4.5 * math.log(2):.6f}   [conjectured exact value of E*(1); b(1) = 3/4]")
    E, m, edges, F = solve_equilibrium(37 / 40, a=6, alpha=3 / 40, n=1000)
    lam, al = 37 / 40, 3 / 40
    Cstar = -2 * lam + 12 * al * lam * (1 - math.log(al)) + 3 * lam ** 2 - 2 * lam ** 2 * math.log(2 * lam)
    print(f"  Fauzan's field (a=6, alpha=3/40, lam=37/40): E* = {E:.5f}, so the true constant U = -E* + C_* = {-E + Cstar:.5f}  (paper proves U <= -1.36700)")


if __name__ == '__main__':
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    part1()
    part3()
    part2([K for K in (10, 20, 30, 40, 50, 60, 80) if K <= Kmax], [(2, 5), (9, 20), (1, 2)])
