#!/usr/bin/env python3
"""
zeta5margin_smodel.py -- Mission B: the asymptotic "local cost model" for the Hankel family
   G = [mu_X^{(s)}(D_N^a t^{i+j}/D_K)]_{i,j<h},  N = alpha K,  h = K - N (lambda = 1 - alpha),  s odd,
against the (numerically exact) equilibrium energy.  Floating point; NOTHING here is a proof.

total_s(a,alpha) = -E*(a,alpha) + C_*(a,alpha) + A_s(a,alpha),   need total < 0 for an irrationality proof.
  E*   : min over measures rho >= 0 of mass lambda on [0,inf) of  int V drho - iint log|y^2-z^2| drho drho,
         V(y) = 2 pi y + int_0^1 log(y^2+u^2)du - a int_0^alpha log(y^2+u^2)du    (s-independent)   [solver of hankel_hadamard.py]
  C_*  : lim K^-2 (log S_K - (K^2 log K)-term),  S_K = (K!)^{2h} 4^{h-1} / ((N!)^{2ah} prod_{i<h}(2i)!^2)
  A_s  : lim K^-2 sum_p (-v_p(S_K) - gamma_p) log p, with the per-prime bound gamma_p of the MODEL:
         x = K/p >= 3 : water filling in the class model, doubled base B = a l_N - l_K - (s-1)
         1 < x < 3    : max( water filling , generalised outer bound (4.12)-(4.14) + Lemma 4.2 )
         x < 1        : 0.
Generalised outer bound (per class with l poles, delta = [class contains a removed pole], needs alpha < 1/6):
   delta=0: cost (s+2)(l-2)^+, zero weights min(l,2);   delta=1: weights min(0, i + a/2 - (l+s-1)/2), i<l-2, and one zero row;
   rank r = (2 lambda - 1 + a alpha) x - (s-1)/2 (per p), loss min(r, #zero weights).
Usage: python3 zeta5margin_smodel.py
"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
import numpy as np
from hankel_hadamard import solve_equilibrium


def types_vec(x, alpha):
    fx = x - np.floor(x); g = alpha * x; gx = g - np.floor(g)
    c = np.stack([fx, 1 - fx, gx, 1 - gx], axis=1)
    c = np.where((c <= 0) | (c >= 0.5), 0.5, c)
    c.sort(axis=1)
    edges = np.concatenate([np.zeros((len(x), 1)), c, np.full((len(x), 1), 0.5)], axis=1)
    mu = np.diff(edges, axis=1)
    z = 0.5 * (edges[:, :-1] + edges[:, 1:])
    X = x[:, None]
    lK = (np.floor(X - z) + np.floor(X + z) + 1).astype(np.int64)
    lN = (np.floor(alpha * X - z) + np.floor(alpha * X + z) + 1).astype(np.int64)
    return mu, lK, lN


def gamma_wf_vec(x, s, a, alpha, lam):
    mu, lK, lN = types_vec(x, alpha)
    B = a * lN - lK - (s - 1)
    rem = lam * x
    tot = np.zeros_like(x)
    w = int(B[mu > 0].min())
    while (rem > 1e-13).any():
        mask = (B <= w) & ((w - B) % 2 == 0)
        take = np.minimum((mu * mask).sum(1), rem)
        tot += take * w
        rem = rem - take
        w += 1
    return tot


def scalN(x, a, alpha, lam):
    m = np.floor(2 * lam * x)
    J = m * lam * x - m * (m + 1) / 4
    return 2 * lam * x * np.floor(x) - 2 * a * lam * x * np.floor(alpha * x) - 2 * J


def cost_out_vec(x, s, a, alpha, lam):
    """generalised (4.14) per p, for 1 < x < 3 and alpha x < 1/2."""
    mu, lK, lN = types_vec(x, alpha)
    cost = np.zeros_like(mu); zeros = np.zeros_like(mu)
    for l in range(0, 8):
        for d in (0, 1):
            sel = (lK == l) & (lN == d)
            if not sel.any(): continue
            if d == 0:
                c, n0 = (s + 2) * max(l - 2, 0), min(l, 2)
            else:
                ws = [min(0.0, i + a / 2 - (l + s - 1) / 2) for i in range(max(l - 2, 0))]
                c = -2 * sum(ws); n0 = sum(1 for w_ in ws if w_ == 0) + (1 if l >= 2 else 0)
            cost[sel] = c; zeros[sel] = n0
    assert (lN[mu > 1e-12] <= 1).all()
    r = np.maximum(0.0, (2 * lam - 1 + a * alpha) * x - (s - 1) / 2)
    return (mu * cost).sum(1) + np.minimum(r, (mu * zeros).sum(1))


def A_model(s, a, alpha, lam, X=40.0, nx=4000, ny=40000, parts=False, outer=True):
    # inner: int_3^X R/x^3, R = -Gamma_wf - N
    x = 3 + (np.arange(int((X - 3) * nx)) + 0.5) / nx
    R = -gamma_wf_vec(x, s, a, alpha, lam) - scalN(x, a, alpha, lam)
    inner = (R / x ** 3).sum() / nx
    # tail: R = x F + Q, F = (s-1) lam + 2 lam {x} - 2 a lam {alpha x};  int_X^inf F/x^2 numerically, Q neglected (|.| <~ 1/X^2)
    xt = X + (np.arange(int(2000 * 400)) + 0.5) / 400
    Ft = (s - 1) * lam + 2 * lam * (xt - np.floor(xt)) - 2 * a * lam * (alpha * xt - np.floor(alpha * xt))
    tail = (Ft / xt ** 2).sum() / 400 + lam * (s - a) / xt[-1]
    # outer: y = p/K in (1/3, 2 lam)
    y = 1 / 3 + (np.arange(ny) + 0.5) * (2 * lam - 1 / 3) / ny
    xo = 1 / y
    C = np.zeros_like(y)
    sel = xo > 1
    cw = -gamma_wf_vec(xo[sel], s, a, alpha, lam)
    if outer and alpha * 3 < 0.5:
        co = cost_out_vec(xo[sel], s, a, alpha, lam)
        C[sel] = np.minimum(cw, co) * y[sel]
        gain = ((co - np.minimum(cw, co)) * y[sel]).sum() * (2 * lam - 1 / 3) / ny
    else:
        C[sel] = cw * y[sel]; gain = float('nan')
    sc = -y * scalN(xo, a, alpha, lam)
    out = (C + sc).sum() * (2 * lam - 1 / 3) / ny
    if parts: return inner, tail, out, gain
    return inner + tail + out


def Cstar(a, alpha, lam):
    t = 2 * a * alpha * lam * (1 - math.log(alpha)) if a * alpha > 0 else 0.0
    return -2 * lam + t + 3 * lam ** 2 - 2 * lam ** 2 * math.log(2 * lam)


def Estar(a, alpha, lam, n=500):
    E, m, edges, F = solve_equilibrium(lam, a=a, alpha=alpha, Y=3.0, n=n)
    sup = edges[1:][m > 1e-12]
    return E, sup.min(), sup.max()


if __name__ == "__main__":
    al, la = 3 / 40, 37 / 40
    i, t, o, g = A_model(5, 6, al, la, parts=True)
    i20 = A_model(5, 6, al, la, X=20.0, parts=True)[0]
    print(f"[check s=5 Fauzan] int_3^20 = {i20:.6f} (exact 0.042788);  int_3^40 + tail = {i + t:.6f} (Fauzan bound for int_3^inf: {0.042788177 - 2689/48000:.6f});"
          f"  outer = {o:.6f} (exact new 1.311113), gain over (4.14) = {g:.6f} (exact 0.019626)")
    E, b0, b = Estar(6, al, la, n=1000)
    print(f"[check s=5 Fauzan] E* = {E:.5f}, C_* = {Cstar(6, al, la):.5f}, U = {-E + Cstar(6, al, la):.5f};  A_5 = {i+t+o:.5f};  total = {-E + Cstar(6, al, la) + i + t + o:+.5f}")
    print("\nscan: total_s(a, alpha) = -E* + C_* + A_s  (lambda = 1 - alpha);  columns s = 3, 5, 7, 9")
    print("   a  alpha |    E*      C_*     U=-E*+C_* [supp]      |   A_3     A_5     A_7     A_9  |  tot_3    tot_5    tot_7    tot_9")
    best = {3: (9, None), 5: (9, None), 7: (9, None), 9: (9, None)}
    grid_a = [0, 2, 4, 5, 6, 7, 8, 9, 10, 12, 14]
    grid_al = [0.02, 0.04, 0.06, 0.075, 0.09, 0.105, 0.12, 0.135, 0.15, 0.165, 0.2, 0.25]
    if len(sys.argv) > 1 and sys.argv[1] == 'fine':
        grid_a = [int(v) for v in sys.argv[2].split(',')]
        grid_al = [float(v) for v in sys.argv[3].split(',')]
    for a in grid_a:
        for alpha in grid_al:
            if a == 0 and alpha != grid_al[0]: continue
            lam = 1 - alpha if a else 1.0
            alp = alpha if a else 1e-9
            E, b0, b = Estar(a, alp, lam)
            cs = Cstar(a, alp, lam)
            A = {s: A_model(s, a, alp, lam, nx=1500, ny=15000) for s in (3, 5, 7, 9)}
            tot = {s: -E + cs + A[s] for s in A}
            for s in A:
                if tot[s] < best[s][0]: best[s] = (tot[s], (a, alpha))
            print(f"  {a:2d} {alpha:6.3f} | {E:8.4f} {cs:8.4f} {-E+cs:9.4f} [{b0:.3f},{b:.3f}] | {A[3]:7.4f} {A[5]:7.4f} {A[7]:7.4f} {A[9]:7.4f} | {tot[3]:+7.4f} {tot[5]:+7.4f} {tot[7]:+8.4f} {tot[9]:+8.4f}" + ("" if alp * 3 < 0.5 else "  (outer: wf only)"), flush=True)
    for s in best:
        print(f"best s={s}: total = {best[s][0]:+.4f} at (a, alpha) = {best[s][1]}")
