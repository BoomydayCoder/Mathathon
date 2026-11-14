#!/usr/bin/env python3
"""
zeta5opt_scan.py -- total(a, alpha) = A(a,alpha) + U(a,alpha) for the family W = D_N^a, N = alpha K, h = K - N, s = 5,
and the resulting irrationality-exponent bound  rho * lambda / eps  (eps = -total).

  A   : zeta5opt_asymptotic.A_numeric  (cost of the proved local bounds, no cutoff; float, error ~1e-6)
  U   : zeta5opt_equilibrium.solve_float  (exact equilibrium energy; float, error ~1e-12)
  rho : Giraudin-type bound for the relative norm ||B_K|| <= e^{rho K}:
        rho = min_{0<a0<b0} [ max_{[a0,b0]} V(t) + max_{alpha<u<=1} ( L(u) + 2 lam arcosh((2u^2+a0+b0)/(b0-a0)) ) ],
        L(u) = a (g(alpha-u) + g(alpha+u)) - g(1-u) - g(1+u),  g(x) = x (log|x| - 1).
        (conditional on Fauzan App. C / Giraudin's Lemmas 3, 5 carrying over to general a -- see the report)
Usage: python3 zeta5opt_scan.py
"""
import sys, math
from fractions import Fraction as Fr
from multiprocessing import Pool
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
import numpy as np
from scipy.optimize import minimize
from zeta5opt_asymptotic import A_numeric
from zeta5opt_equilibrium import solve_float, Vt


def rho_norm(a, alpha, lam, start=(1e-3, 0.245)):
    g = lambda x: np.where(x == 0, 0.0, x * (np.log(np.abs(np.where(x == 0, 1, x))) - 1))
    u = alpha + (1 - alpha) * (np.arange(1, 6001) / 6000)
    L = a * (g(alpha - u) + g(alpha + u)) - g(1 - u) - g(1 + u)

    def val(z):
        a0, b0 = math.exp(z[0]), math.exp(z[0]) + math.exp(z[1])
        t = np.exp(np.linspace(math.log(a0), math.log(b0), 1500))
        Phi = L + 2 * lam * np.arccosh((2 * u * u + a0 + b0) / (b0 - a0))
        return Vt(t, a, alpha).max() + Phi.max()
    r = minimize(val, [math.log(start[0]), math.log(start[1] - start[0])], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-10, maxiter=4000))
    a0, b0 = math.exp(r.x[0]), math.exp(r.x[0]) + math.exp(r.x[1])
    return r.fun, a0, b0, val([math.log(start[0]), math.log(start[1] - start[0])])


def job(arg):
    a, num, den = arg
    alpha, lam = num / den, 1 - num / den
    r = A_numeric(a, num, den, X0=40, nx=2000, ny=30000)
    e = solve_float(a, alpha, lam)
    rho, a0, b0, rho_g = rho_norm(a, alpha, lam)
    tot = r['A'] + e['U']
    return (a, num, den, r['A'], e['U'], tot, e['A'], e['B'], rho, a0, b0, rho * lam / (-tot) if tot < 0 else float('inf'), r['gain_wf'])


if __name__ == "__main__":
    rg, _, _, rG = rho_norm(6, 3 / 40, 37 / 40)
    print(f"[check] rho at Fauzan's parameters with Giraudin's (a0,b0) = (1e-3, 49/200): {rG:.7f}  (Giraudin: 9.0827901);  optimised over (a0,b0): {rg:.7f}")
    alphas = [(1, 16), (1, 14), (3, 40), (1, 12), (1, 11), (1, 10), (2, 19), (1, 9), (9, 80), (2, 17), (3, 25), (1, 8), (2, 15), (1, 7), (3, 20)]
    jobs = [(a, n, d) for a in (4, 5, 6, 7, 8) for n, d in alphas]
    with Pool(4) as P:
        res = P.map(job, jobs)
    print("  a  alpha        |    A         U        total    | supp sqrt[A,B]     | rho      (a0, b0)          | exponent rho*lam/eps | wf gain in outer")
    for (a, n, d, A, U, tot, sA, sB, rho, a0, b0, ex, gw) in res:
        print(f"  {a}  {n:2d}/{d:<3d}={n/d:.4f} | {A:9.6f} {U:10.6f} {tot:+9.6f} | [{math.sqrt(sA):.4f},{math.sqrt(sB):.4f}] | {rho:7.4f} ({a0:.5f},{b0:.4f}) | {ex:8.2f} | {gw:.5f}")
    best = min(res, key=lambda r: r[5])
    bex = min(res, key=lambda r: r[11])
    print(f"best total: {best[5]:+.6f} at a={best[0]}, alpha={best[1]}/{best[2]}")
    print(f"best exponent: {bex[11]:.2f} at a={bex[0]}, alpha={bex[1]}/{bex[2]}  (total {bex[5]:+.6f}, rho {bex[8]:.4f})")
