#!/usr/bin/env python3
"""Verifier's independent trace-LP (own constraint generation, own enumeration) on FRESH Goldstein-Mayer lattices with a
LARGE modulus p (non-granular), reduced by a QR-guided integer LLL.  Usage: fresh_lp.py n p nsamp timecap seed out"""
import sys, os, json, math, time
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import *

def int_lll(B, delta=0.99):
    """LLL on the columns of an integer matrix (int64), guided by a float QR recomputed after every swap."""
    B = B.copy(); n = B.shape[1]
    R = np.linalg.qr(B.astype(float), mode='r'); k = 1; it = 0
    while k < n:
        it += 1
        for j in range(k - 1, -1, -1):
            q = int(np.rint(R[j, k] / R[j, j]))
            if q:
                B[:, k] -= q * B[:, j]; R[:, k] -= q * R[:, j]
        if R[k, k] ** 2 + R[k - 1, k] ** 2 >= delta * R[k - 1, k - 1] ** 2:
            k += 1
        else:
            B[:, [k - 1, k]] = B[:, [k, k - 1]]
            R = np.linalg.qr(B.astype(float), mode='r'); k = max(k - 1, 1)
    assert np.abs(B).max() < 2 ** 40
    return B

def gm_lattice(n, p, rng):
    B = np.eye(n, dtype=np.int64)
    B[n - 1, :n - 1] = rng.integers(0, p, size=n - 1); B[n - 1, n - 1] = p
    B = int_lll(B); B = int_lll(B)
    G = (B.T @ B)                                   # exact integer Gram matrix
    d = round(abs(np.linalg.det(B.astype(float))) / p, 6)
    assert d == 1.0, d                              # index p sublattice of Z^n
    Q = G / p ** (2.0 / n)
    Q, U = my_lll(Q)
    return Q / np.linalg.det(Q) ** (1.0 / n), int(G.diagonal().min())

def canon(X):
    first = np.argmax(X != 0, axis=1); sg = np.sign(X[np.arange(len(X)), first]); return X * sg[:, None]

def my_trace_lp(Q):
    n = Q.shape[0]; N = n * (n + 1) // 2; iu = np.triu_indices(n); w = np.where(iu[0] == iu[1], 1.0, 2.0)
    Qi = np.linalg.inv(Q); c = Qi[iu] * w
    X, v, _ = short(Q, Q.diagonal().min() * (1 + 1e-12)); lam2 = float(v.min())
    bound = 1.4 * lam2
    while True:
        X, v, _ = short(Q, bound, maxc=400000)
        if len(X) >= 6 * N: break
        bound *= 1.12
    cons = {tuple(x) for x in canon(X).tolist()}
    rounds = 0
    while True:
        rounds += 1
        Xc = np.array(sorted(cons), dtype=np.int64); Xf = Xc.astype(float)
        Phi = Xf[:, iu[0]] * Xf[:, iu[1]] * w
        res = linprog(c, A_ub=-Phi, b_ub=-np.ones(len(Xc)), bounds=[(-1e4, 1e4)] * N, method='highs')
        assert res.status == 0, res.message
        A = np.zeros((n, n)); A[iu] = res.x; A = A + A.T - np.diag(np.diag(A))
        ev = np.linalg.eigvals(Qi @ A).real
        bad = np.abs(res.x).max() > 0.99e4 or ev.min() < 1e-4 * ev.max()
        if not bad:
            try:
                Xv, vv, _ = short(A, 1 - 1e-9, maxc=400000, maxnodes=20_000_000)
            except RuntimeError:
                bad = True
        if bad:
            bound *= 1.2; X, v, _ = short(Q, bound, maxc=400000); cons |= {tuple(x) for x in canon(X).tolist()}; continue
        if len(Xv) == 0: break
        new = [tuple(x) for x in canon(Xv[np.argsort(vv)]).tolist() if tuple(x) not in cons][:4 * N]
        assert new, 'violators already in the constraint set (LP tolerance)'
        cons |= set(new)
        assert rounds < 500
    S, vs, _ = short(A, 1 + 1e-7)
    M = np.einsum('ki,kj->kij', S, S).reshape(len(S), n * n).astype(float)
    dl = np.linalg.lstsq(M, 1.0 - M @ A.ravel(), rcond=None)[0].reshape(n, n); A = A + (dl + dl.T) / 2
    return A, lam2, rounds, len(cons)

if __name__ == '__main__':
    n = int(sys.argv[1]); p = int(sys.argv[2]); nsamp = int(sys.argv[3]); cap = float(sys.argv[4]); seed = int(sys.argv[5]); out = sys.argv[6]
    rng = np.random.default_rng([n, seed, 424242]); N = n * (n + 1) // 2; recs = []; t0 = time.time()
    for i in range(nsamp):
        if time.time() - t0 > cap: break
        Q, gmin = gm_lattice(n, p, rng)
        A, lam2, rounds, ncons = my_trace_lp(Q)
        mn, npairs, S = min_and_pairs(A)
        V = kap(n) / math.sqrt(np.linalg.det(A)); Qi = np.linalg.inv(Q)
        rr, sl, npos, lam = dual_certificate(Q, S)
        iu = np.triu_indices(n); Phi = (S[:, iu[0]] * S[:, iu[1]]).astype(float)
        resid = float(np.linalg.norm(Phi.T @ lam - Qi[iu]) / np.linalg.norm(Qi[iu]))
        m = np.trace(Qi @ A) / n
        assert abs(mn - 1) < 1e-9 and V <= Vbest(n) * (1 + 1e-7) and resid < 1e-8, (mn, V, resid)
        recs.append(dict(n=n, p=p, idx=i, grid_lam1sq=gmin, lam1sq=lam2, V0=kap(n) * lam2 ** (n / 2), V_LP=V, v_LP=kap(n) * m ** (-n / 2), pairs=npairs, min=mn,
                         dual_resid=resid, rounds=rounds, ncons=ncons, best=bool(V >= Vbest(n) * (1 - 1e-6)), t=time.time() - t0))
        json.dump(recs, open(out, 'w'))
    V = np.array([r['V_LP'] for r in recs]) / N
    print(f'fresh n={n} p={p}: {len(recs)} lattices in {time.time()-t0:.1f}s | V_LP/N = {V.mean():.3f}+-{V.std(ddof=1)/math.sqrt(len(V)):.3f} (std {V.std(ddof=1):.3f}) | v_LP/N = {np.mean([r["v_LP"] for r in recs])/N:.3f} | '
          f'pairs/N = {np.mean([r["pairs"] for r in recs])/N:.3f} | simple {sum(r["pairs"]==N for r in recs)}/{len(recs)} | best known {sum(r["best"] for r in recs)}/{len(recs)} | mean V0 = {np.mean([r["V0"] for r in recs]):.2f} | '
          f'mean grid lam1^2 = {np.mean([r["grid_lam1sq"] for r in recs]):.1f} | max dual resid {max(r["dual_resid"] for r in recs):.1e} | max |min-1| {max(abs(r["min"]-1) for r in recs):.1e}', flush=True)
