#!/usr/bin/env python3
"""Exact trace-LP for a Poisson process of pairs in R^n (no truncation error): row generation.
Points x = t^{1/n} u (units: ball through x has volume t), pairs have rate 1/2 per unit volume.
LP: minimise Tr A s.t. x^T A x >= 1 for all points.  Start with points of volume-radius t <= T; then repeatedly
sample the Poisson process in (E_k minus B_T minus earlier E_j) (independent of everything sampled before),
add the points found as constraints and re-solve, until the current ellipsoid is empty.
The final A is optimal for the whole (infinite) process.  v = (Tr A/n)^{-n/2}, V = det(A)^{-1/2}."""
import sys, time
import numpy as np
from scipy.optimize import linprog

def run(n, seed, srel=9.0, maxit=60):
    rng = np.random.default_rng(seed)
    N = n * (n + 1) // 2
    T = (2 * N / np.e) * np.exp(srel / 2)
    M = rng.poisson(T / 2)
    t = T * rng.random(M)
    U = rng.standard_normal((M, n)); U /= np.linalg.norm(U, axis=1, keepdims=True)
    iu = np.triu_indices(n); w = np.where(iu[0] == iu[1], 1.0, 2.0)
    F = U[:, iu[0]] * U[:, iu[1]] * w
    rhs = t ** (-2.0 / n)
    c = np.where(iu[0] == iu[1], 1.0, 0.0)
    prev = []; added = 0; t0 = time.time()
    for it in range(maxit):
        res = linprog(c, A_ub=-F, b_ub=-rhs, bounds=[(None, None)] * N, method="highs")
        if res.status != 0:
            return dict(n=n, seed=seed, status=res.status)
        Am = np.zeros((n, n)); Am[iu] = res.x; Am = Am + Am.T - np.diag(np.diag(Am))
        ev, Q = np.linalg.eigh(Am)
        if ev[0] <= 0:
            return dict(n=n, seed=seed, status="not PD", it=it)
        V = float(np.prod(ev) ** (-0.5))
        K = rng.poisson(V / 2)
        Y = rng.standard_normal((K, n)); Y /= np.linalg.norm(Y, axis=1, keepdims=True)
        Y *= rng.random((K, 1)) ** (1.0 / n)
        X = (Y @ (Q * ev ** (-0.5)).T)            # x = A^{-1/2} y, uniform in E_k
        r = np.linalg.norm(X, axis=1)
        keep = r ** n > T
        for Ap in prev:
            keep &= np.einsum("ij,jk,ik->i", X, Ap, X) >= 1.0
        prev.append(Am)
        Xn = X[keep]
        if len(Xn) == 0:
            break
        Un = Xn / np.linalg.norm(Xn, axis=1, keepdims=True); tn = np.linalg.norm(Xn, axis=1) ** n
        F = np.vstack([F, Un[:, iu[0]] * Un[:, iu[1]] * w]); rhs = np.concatenate([rhs, tn ** (-2.0 / n)])
        added += len(Xn)
    else:
        return dict(n=n, seed=seed, status="no convergence", added=added)
    tr = np.trace(Am) / n
    v = tr ** (-n / 2.0)
    B0 = Am / tr - np.eye(n); rho2 = 2 * float(np.sum(B0 * B0))
    nact = int(np.sum(F @ res.x - rhs < 1e-9))
    return dict(n=n, N=N, seed=seed, iters=it + 1, added=added, time=round(time.time() - t0, 1), v_over_N=round(v / N, 4),
                rho2=round(rho2, 3), V_over_N=round(V / N, 4), V_over_n2=round(V / n**2, 4), V_2nd_over_N=round(v * np.exp(rho2 / 8) / N, 4),
                act_over_N=nact / N, lam_min=round(ev[0] / tr, 4), lam_max=round(ev[-1] / tr, 4))

if __name__ == "__main__":
    for j in sys.argv[1:]:
        n, seeds = j.split(":")
        for s in range(int(seeds)):
            print(run(int(n), 3000 + s), flush=True)
