#!/usr/bin/env python3
"""Numerical check of claim (C1) in the Gaussian/Poisson MODEL and in the exact Poisson-Veronese model.

Model G (Gaussian): constraints <theta,g_i> >= -sigma_i, g_i iid N(0,I_d), d=N-1,
  sigma_i a Poisson process with intensity (v/4) e^{sigma/2}.  Writing sigma = sigma' - 2 log v with
  sigma' of intensity (1/4) e^{sigma'/2}, "largest v with K_v nonempty" is the LP
      maximize tau  s.t.  <theta,g_i> - tau >= -sigma'_i ,   v* = exp(tau*/2).
  Prediction (C1): v*/N -> 2/e = 0.7358, |theta*|^2 -> 8, V = v* e^{|theta*|^2/8} -> 2N, #active = N.

Model P (exact Poisson points in R^n, no Gaussian approximation, no linearisation):
  pairs +-x_i, x_i = t_i^{1/n} u_i (units: ball through x_i has volume t_i), t_i Poisson process of rate 1/2
  on (0,T), u_i uniform on the sphere.  LP: minimise Tr A s.t. x_i^T A x_i >= 1.
  v = (Tr A/n)^{-n/2} (trace-equivalent ball), V = det(A)^{-1/2} (true volume of the empty ellipsoid).
  A posteriori exactness of the truncation: the ellipsoid is contained in the ball of volume
  lambda_min(A)^{-n/2}; if that is <= T no point beyond T can lie inside.
"""
import sys, time
import numpy as np
from scipy.optimize import linprog

def run_gauss(d, seed, srel=11.0, extra=3.0):
    rng = np.random.default_rng(seed)
    N = d + 1
    vpred = 2 * N / np.e
    S = srel + 2 * np.log(vpred)
    M = rng.poisson(0.5 * np.exp(S / 2))
    sp = S + 2 * np.log(rng.random(M))
    G = rng.standard_normal((M, d))
    A = np.hstack([-G, np.ones((M, 1))])
    c = np.zeros(d + 1); c[-1] = -1.0
    t0 = time.time()
    res = linprog(c, A_ub=A, b_ub=sp, bounds=[(None, None)] * (d + 1), method="highs")
    tt = time.time() - t0
    if res.status != 0:
        return dict(model="G", d=d, seed=seed, status=res.status, msg=res.message)
    th = res.x[:d]; tau = res.x[-1]
    v = np.exp(tau / 2); rho2 = float(th @ th)
    slack = sp - A @ res.x
    nact = int(np.sum(slack < 1e-7))
    # check constraints beyond the truncation: sigma' in (S, S+extra)
    M2 = rng.poisson(0.5 * (np.exp((S + extra) / 2) - np.exp(S / 2)))
    viol = 0; done = 0
    while done < M2:
        m = min(20000, M2 - done)
        u = rng.random(m)
        s2 = 2 * np.log(np.exp(S / 2) + u * (np.exp((S + extra) / 2) - np.exp(S / 2)))
        g2 = rng.standard_normal((m, d))
        viol += int(np.sum(g2 @ th - tau < -s2))
        done += m
    return dict(model="G", d=d, N=N, seed=seed, M=M, time=round(tt, 1), v_over_N=v / N, rho2=rho2,
                V_over_N=v * np.exp(rho2 / 8) / N, act_over_N=nact / N, viol_beyond=viol, M_beyond=M2)

def run_poisson(n, seed, srel=11.0):
    rng = np.random.default_rng(seed)
    N = n * (n + 1) // 2
    vpred = 2 * N / np.e
    T = vpred * np.exp(srel / 2)
    M = rng.poisson(T / 2)
    t = T * rng.random(M)
    U = rng.standard_normal((M, n)); U /= np.linalg.norm(U, axis=1, keepdims=True)
    iu = np.triu_indices(n)
    w = np.where(iu[0] == iu[1], 1.0, 2.0)          # x^T A x = sum_{i<=j} w_ij a_ij x_i x_j
    F = U[:, iu[0]] * U[:, iu[1]] * w                 # M x N
    rhs = t ** (-2.0 / n)                             # u^T A u >= t^{-2/n}
    c = np.where(iu[0] == iu[1], 1.0, 0.0)            # Tr A
    t0 = time.time()
    res = linprog(c, A_ub=-F, b_ub=-rhs, bounds=[(None, None)] * N, method="highs")
    tt = time.time() - t0
    if res.status != 0:
        return dict(model="P", n=n, seed=seed, status=res.status, msg=res.message)
    Am = np.zeros((n, n)); Am[iu] = res.x; Am = Am + Am.T - np.diag(np.diag(Am))
    ev = np.linalg.eigvalsh(Am)
    v = (np.trace(Am) / n) ** (-n / 2.0)
    V = float(np.prod(ev) ** (-0.5)) if ev[0] > 0 else float("nan")
    B0 = Am / (np.trace(Am) / n) - np.eye(n)
    rho2 = 2 * float(np.sum(B0 * B0))
    slack = F @ res.x - rhs
    nact = int(np.sum(slack < 1e-9))
    circ = ev[0] ** (-n / 2.0) if ev[0] > 0 else float("inf")   # volume of circumscribed ball
    return dict(model="P", n=n, N=N, seed=seed, M=M, time=round(tt, 1), v_over_N=v / N, rho2=rho2,
                V_over_N=V / N, V_2nd_order_over_N=v * np.exp(rho2 / 8) / N, act_over_N=nact / N,
                lam_min=ev[0] / (np.trace(Am) / n), lam_max=ev[-1] / (np.trace(Am) / n),
                trunc_exact=bool(circ <= T), circ_over_T=circ / T)

if __name__ == "__main__":
    jobs = sys.argv[1:]
    for j in jobs:
        kind, size, seeds = j.split(":")
        for s in range(int(seeds)):
            r = run_gauss(int(size), 1000 + s) if kind == "G" else run_poisson(int(size), 2000 + s)
            print({k: (round(x, 4) if isinstance(x, float) else x) for k, x in r.items()}, flush=True)
