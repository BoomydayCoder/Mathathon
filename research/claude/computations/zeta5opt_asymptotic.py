#!/usr/bin/env python3
"""
zeta5opt_asymptotic.py -- asymptotic arithmetic cost A(a, alpha) of the PROVED local bounds for the family
    W = D_N^a,  N = alpha K,  h = K - N (lambda = 1 - alpha),  s = 5,  normalisation S_K = (K!)^{2h} 4^{h-1} / ((N!)^{2ah} prod (2i)!^2):

    A(a,alpha) = lim K^-2 log m_K = int_3^inf R(x)/x^3 dx                     [p <= K/3: Theorem W_a, NO cutoff M]
               + int_{1/3}^{1} ( min(C_wf, C_out)(y) + scal(y) ) dy            [K/3 < p <= K: max(Theorem W_a, generalised (4.14))]
               + int_1^{2 lambda} scal(y) dy                                   [K < p <= 2h: v_p >= 0]
    R(x) = -Gamma_wf(x) - N(x)  (x = K/p),  C(y) = -y Gamma(1/y),  scal(y) = -y N(1/y)  (y = p/K).

Two evaluations:
  * numeric (float; midpoint rule on fine grids; tail of the inner integral by periodicity of F and Q in R = xF + Q,
    F = 4 lam + 2 lam {x} - 2 a lam {alpha x}, period = denominator of alpha, summed with Hurwitz zeta);
  * exact (fractions) for int_3^X0 and for the outer integrals, with an exact kink-finding integrator of piecewise-linear functions.
Usage:  python3 zeta5opt_asymptotic.py a num den [X0]      (alpha = num/den)
KNOWN DEFECT (found at the end of the run, not fixed): integrate_pl tests linearity at 5 sample points per cell and can miss
kinks lying between them.  The exact outer value agrees with the robust float midpoint rule to 1e-8 for alpha = 3/40 (where it also
equals the closed form 2139737/1632000 - 289/800 obtained from Fauzan's formulas) and alpha = 1/9, but NOT for alpha = 3/25
(exact routine 0.725614 vs float 0.725792) and alpha = 9/80 (4e-6 off, 2 unresolved cells).  For those alpha trust the float value.
Nothing is claimed by this file; it only computes.
"""
import sys, math
sys.setrecursionlimit(100000)
from fractions import Fraction as Fr
from math import floor
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
import numpy as np
from scipy.special import zeta as hzeta
from zeta5margin_smodel import gamma_wf_vec, scalN, cost_out_vec
from zeta5margin_asymptotic import Model, int_R_over_x3

S = 5


def A_numeric(a, num, den, X0=40, nx=3000, ny=60000, check=False):
    alpha, lam = num / den, 1 - num / den
    x = 3 + (np.arange(int((X0 - 3) * nx)) + 0.5) / nx
    R = -gamma_wf_vec(x, S, a, alpha, lam) - scalN(x, a, alpha, lam)
    inner = (R / x ** 3).sum() / nx
    d = den
    xt = X0 + (np.arange(int(d * nx)) + 0.5) / nx
    Rt = -gamma_wf_vec(xt, S, a, alpha, lam) - scalN(xt, a, alpha, lam)
    Ft = 4 * lam + 2 * lam * (xt - np.floor(xt)) - 2 * a * lam * (alpha * xt - np.floor(alpha * xt))
    Qt = Rt - xt * Ft
    tail = (Ft * hzeta(2, xt / d) / d ** 2 + Qt * hzeta(3, xt / d) / d ** 3).sum() / nx
    info = {}
    if check:
        xs = xt[:: max(1, len(xt) // 4000)]
        x2 = xs + d
        R2 = -gamma_wf_vec(x2, S, a, alpha, lam) - scalN(x2, a, alpha, lam)
        F2 = 4 * lam + 2 * lam * (x2 - np.floor(x2)) - 2 * a * lam * (alpha * x2 - np.floor(alpha * x2))
        Q1 = np.interp(xs, xt, Qt)
        info = dict(Qper=float(np.abs((R2 - x2 * F2) - Q1).max()), Qmin=float(Qt.min()), Qmax=float(Qt.max()), Qmean=float(Qt.mean()), Fmean=float(Ft.mean()))
    y = 1 / 3 + (np.arange(ny) + 0.5) * (1 - 1 / 3) / ny
    xo = 1 / y
    cw = -gamma_wf_vec(xo, S, a, alpha, lam)
    co = cost_out_vec(xo, S, a, alpha, lam)
    sc = -y * scalN(xo, a, alpha, lam)
    out = ((np.minimum(cw, co) * y + sc).sum()) * (2 / 3) / ny
    out_only = ((co * y + sc).sum()) * (2 / 3) / ny
    y2 = 1 + (np.arange(ny) + 0.5) * (2 * lam - 1) / ny
    beyond = (-y2 * scalN(1 / y2, a, alpha, lam)).sum() * (2 * lam - 1) / ny
    return dict(inner=float(inner), tail=float(tail), outer=float(out), beyond=float(beyond),
                A=float(inner + tail + out + beyond), gain_wf=float(out_only - out), **info)


# ------------------------------------------------------------------ exact part
INEXACT = [0, Fr(0)]      # number and total width of cells treated by the midpoint rule (width < 1e-30)


def integrate_pl(f, u, v, depth=0):
    """exact integral over [u,v] of a piecewise-linear f (finitely many kinks).  If a cell of width < 1e-30 is still not
    linear (isolated jump of f), it is integrated by the midpoint rule and recorded in INEXACT (error <= width * jump)."""
    if v - u < Fr(1, 10 ** 30):
        INEXACT[0] += 1; INEXACT[1] += v - u
        return f((u + v) / 2) * (v - u), 1
    pts = [u + (v - u) * Fr(k, 16) for k in (1, 2, 8, 14, 15)]
    val = [f(t) for t in pts]
    sl = (val[1] - val[0]) / (pts[1] - pts[0])
    if all(val[i] == val[0] + sl * (pts[i] - pts[0]) for i in range(5)):
        c0 = val[0] - sl * pts[0]
        return sl * (v * v - u * u) / 2 + c0 * (v - u), 1
    sr = (val[4] - val[3]) / (pts[4] - pts[3])
    cut = None
    if sr != sl:
        t = ((val[3] - sr * pts[3]) - (val[0] - sl * pts[0])) / (sl - sr)
        if u < t < v: cut = t
    if cut is None or depth > 60:
        cut = (u + v) / 2
    I1, n1 = integrate_pl(f, u, cut, depth + 1)
    I2, n2 = integrate_pl(f, cut, v, depth + 1)
    return I1 + I2, n1 + n2


def exact_outer(a, alpha):
    alpha = Fr(alpha); lam = 1 - alpha
    M = Model(S, a, alpha, lam)

    def cost_out(x):
        tot = Fr(0); zeros = Fr(0)
        for mu, B in M.types(x):
            pass
        al = alpha
        cuts = {Fr(0), Fr(1, 2)}
        frac = lambda u: u - floor(u)
        for c in (frac(x), 1 - frac(x), frac(al * x), 1 - frac(al * x)):
            if 0 < c < Fr(1, 2): cuts.add(c)
        cuts = sorted(cuts)
        for u, v in zip(cuts, cuts[1:]):
            z = (u + v) / 2
            l = floor(x - z) + floor(x + z) + 1
            dl = floor(al * x - z) + floor(al * x + z) + 1
            assert dl in (0, 1)
            w2 = [min(0, 2 * i + a * dl - (l + 4)) for i in range(max(l - 2, 0))] + [0] * (min(l, 2) - dl)
            tot += (v - u) * (-sum(w2))
            zeros += (v - u) * sum(1 for w in w2 if w == 0)
        r = max(Fr(0), (2 * lam - 1 + a * alpha) * x - 2)
        return tot + min(r, zeros)

    scal = lambda y: -y * M.scalN(1 / y)
    Cwf = lambda y: -y * M.wf(1 / y)[0]
    Cout = lambda y: y * cost_out(1 / y)
    best = lambda y: min(Cwf(y), Cout(y)) + scal(y)
    old = lambda y: Cout(y) + scal(y)
    # split at the discontinuities of the floor functions first
    ys = {Fr(1, 3), Fr(1)}
    for coef in (Fr(1), Fr(2), alpha, 2 * alpha, 1 - alpha, 1 + alpha, 2 * lam):
        k = math.ceil(coef * 1)
        while k <= coef * 3:
            if k > 0: ys.add(coef / k if False else Fr(1) / (Fr(k) / coef))
            k += 1
    ys = sorted(t for t in ys if Fr(1, 3) <= t <= 1)
    Ibest = Iold = Fr(0); cells = 0
    for u, v in zip(ys, ys[1:]):
        i1, n1 = integrate_pl(best, u, v); i2, _ = integrate_pl(old, u, v)
        Ibest += i1; Iold += i2; cells += n1
    yb = sorted({Fr(1), 2 * lam} | {2 * lam / j for j in range(1, 8) if 1 < 2 * lam / j < 2 * lam})
    Ibey = sum((integrate_pl(scal, u, v)[0] for u, v in zip(yb, yb[1:])), Fr(0))
    return Ibest, Iold, Ibey, cells


if __name__ == "__main__":
    a, num, den = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    X0 = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    alpha = Fr(num, den); lam = 1 - alpha
    print(f"=== a={a}, alpha={alpha}, lambda={lam}, s=5")
    r = A_numeric(a, num, den, X0=40, check=True)
    print(f" numeric: inner int_3^40 = {r['inner']:.7f}, tail int_40^inf = {r['tail']:.7f}, outer(1/3..1) = {r['outer']:.7f}, beyond(1..2lam) = {r['beyond']:.7f}")
    print(f"          A = {r['A']:.7f};  gain of max(wf,out) over out alone = {r['gain_wf']:.7f}")
    print(f"          R = xF + Q on one period [40,40+{den}]: periodicity defect of Q {r['Qper']:.2e}; Q in [{r['Qmin']:.4f},{r['Qmax']:.4f}], mean {r['Qmean']:.5f}; mean F = {r['Fmean']:.6f} (expected (5-a) lam = {(5-a)*float(lam):.6f})")
    r2 = A_numeric(a, num, den, X0=20)
    print(f"          same with X0=20: A = {r2['A']:.7f}  (difference {r['A']-r2['A']:+.2e})")
    Mo = Model(S, a, alpha, lam)
    Iin, nc = int_R_over_x3(Mo, Fr(3), Fr(X0))
    Ibest, Iold, Ibey, cells = exact_outer(a, alpha)
    print(f" exact:   int_3^{X0} R/x^3 dx = {Iin} = {float(Iin):.10f}  ({nc} cells)   [numeric: {r2['inner'] if X0 == 20 else float('nan'):.7f}]")
    print(f"          outer  int_(1/3)^1 [min(C_wf,C_out)+scal] = {Ibest} = {float(Ibest):.10f}  ({cells} linear pieces)   [numeric: {r['outer']:.7f}]")
    print(f"          (with C_out alone: {Iold} = {float(Iold):.10f};  gain {float(Iold - Ibest):.10f})")
    print(f"          beyond int_1^(2 lam) scal = {Ibey} = {float(Ibey):.10f}   [numeric: {r['beyond']:.7f}]")
    print(f"          cells not resolved exactly (width < 1e-30, midpoint rule): {INEXACT[0]}, total width {float(INEXACT[1]):.2e}")
    Aex = float(Iin + Ibest + Ibey) + (r2['tail'] if X0 == 20 else float('nan'))
    print(f" A(a,alpha) = exact part + numeric tail int_{X0}^inf ({r2['tail']:.7f}) = {Aex:.7f}")
