#!/usr/bin/env python3
"""Simulation of Klartag's stochastically evolving ellipsoid (sticky contact points).
dA = pi_F(dW), W Brownian motion on Sym_n with <A,B>=Tr(AB); F = {B: Bx.x=0 for contact x}.
Usage: python3 sim.py n kind(gm|zn) maxruns timecap_seconds seed outfile
"""
import sys, json, time, math, os
import numpy as np
from math import gamma, pi

def kappa(n): return pi ** (n / 2) / gamma(n / 2 + 1)

def lll(M, extras=(), delta=0.99):
    """LLL on columns of M (in place, column ops also applied to matrices in extras)."""
    n = M.shape[1]
    k = 1
    while k < n:
        G = M.T @ M
        R = np.linalg.cholesky(G).T
        for j in range(k - 1, -1, -1):
            q = np.rint(R[j, k] / R[j, j])
            if q != 0:
                R[:, k] -= q * R[:, j]
                M[:, k] -= q * M[:, j]
                for E in extras: E[:, k] -= q * E[:, j]
        if R[k, k] ** 2 + R[k - 1, k] ** 2 >= delta * R[k - 1, k - 1] ** 2:
            k += 1
        else:
            M[:, [k - 1, k]] = M[:, [k, k - 1]]
            for E in extras: E[:, [k - 1, k]] = E[:, [k, k - 1]]
            k = max(k - 1, 1)
    return M

def enum_points(G, bound):
    """All nonzero z in Z^n with z^T G z <= bound, one representative per +-pair (Fincke-Pohst, BFS)."""
    n = G.shape[0]
    R = np.linalg.cholesky(G).T
    Z = np.zeros((1, n)); P = np.zeros(1)
    b = bound * (1 + 1e-12)
    for i in range(n - 1, -1, -1):
        c = -(Z[:, i + 1:] @ R[i, i + 1:]) / R[i, i]
        rad = np.sqrt(np.maximum(b - P, 0.0)) / R[i, i]
        lo = np.ceil(c - rad); hi = np.floor(c + rad)
        if i == n - 1: lo = np.maximum(lo, 0)
        cnt = (hi - lo + 1).astype(np.int64); cnt[cnt < 0] = 0
        tot = int(cnt.sum())
        idx = np.repeat(np.arange(len(cnt)), cnt)
        starts = np.cumsum(cnt) - cnt
        off = np.arange(tot) - np.repeat(starts, cnt)
        zi = lo[idx] + off
        ci = c[idx]
        Z = Z[idx]; Z[:, i] = zi
        P = P[idx] + (R[i, i] * (zi - ci)) ** 2
    nz = Z != 0
    anynz = nz.any(1)
    last = n - 1 - np.argmax(nz[:, ::-1], axis=1)
    sg = Z[np.arange(len(Z)), last]
    keep = anynz & (sg > 0)
    return Z[keep]

def quad(A, X): return ((X @ A) * X).sum(1)

class Sim:
    def __init__(self, n, B, rng, dt, eta=1.0):
        self.n = n; self.d = n * (n + 1) // 2; self.rng = rng; self.dt = dt
        self.eta = eta; self.thr = eta / (1 + eta)
        self.iu = np.triu_indices(n, 1); self.di = np.diag_indices(n)
        self.covol = abs(np.linalg.det(B)); self.kap = kappa(n)
        self.B = B.copy()
        lll(self.B)
        G = self.B.T @ self.B
        bd = G.diagonal().min()
        Z = enum_points(G, bd)
        self.lam1 = math.sqrt(quad(G, Z).min())   # lambda_1 in the covol = kappa_n normalization
        # scale invariance (A -> cA  <=>  t -> c^2 t): rescale the lattice so that lambda_1 = 1 and a_0 = 1.001,
        # i.e. Klartag's normalization a_0 ~ 1; V is computed with the actual covolume, so it is unchanged.
        self.B = self.B / self.lam1
        self.covol = abs(np.linalg.det(self.B))
        self.a0 = 1.001
        self.A = self.a0 * np.eye(n)
        self.contacts = np.zeros((0, n)); self.Qc = np.zeros((self.d, 0))
        self.t = 0.0; self.nenum = 0; self.warn = 0
        self.events = []  # (t, V, npairs, rank)
        self.reenum()

    def vec(self, X):
        return np.concatenate([X * X, math.sqrt(2) * X[:, self.iu[0]] * X[:, self.iu[1]]], axis=1)
    def vecA(self, A):
        return np.concatenate([A[self.di], math.sqrt(2) * A[self.iu]])
    def unvec(self, v):
        n = self.n
        D = np.zeros((n, n)); D[self.iu] = v[n:] / math.sqrt(2); D = D + D.T; D[self.di] = v[:n]
        return D
    def V(self):
        return self.kap / math.sqrt(np.linalg.det(self.A)) / self.covol
    def add_contacts(self, Xn):
        self.contacts = np.vstack([self.contacts, Xn])
        for v in self.vec(Xn):
            w = v.copy()
            for _ in range(2):
                w -= self.Qc @ (self.Qc.T @ w)
            nr = np.linalg.norm(w)
            if nr > 1e-8 * np.linalg.norm(v):
                self.Qc = np.hstack([self.Qc, (w / nr)[:, None]])
    def fix_drift(self):
        if len(self.contacts) == 0: return
        Mc = self.vec(self.contacts)
        r = 1.0 - Mc @ self.vecA(self.A)
        vS = np.linalg.lstsq(Mc, r, rcond=None)[0]
        self.A = self.A + self.unvec(vS)
    def reenum(self):
        self.nenum += 1
        self.fix_drift()
        Lc = np.linalg.cholesky(self.A)
        M = Lc.T @ self.B
        lll(M, extras=(self.B,))
        Z = enum_points(M.T @ M, 1 + self.eta)
        X = Z @ self.B.T
        q = quad(self.A, X)
        if len(q) and q.min() < 1 - 1e-7: self.warn += 1
        near = q <= 1 + 1e-9
        if near.sum() != len(self.contacts):
            if near.sum() < len(self.contacts): self.warn += 1
            # rebuild contact set from enumeration
            self.contacts = np.zeros((0, self.n)); self.Qc = np.zeros((self.d, 0))
            self.add_contacts(X[near])
        self.X = X[~near]; self.q = q[~near]
        self.A0 = self.A.copy(); self.L0inv = np.linalg.inv(Lc)
    def rel(self, S):
        return np.linalg.norm(self.L0inv @ S @ self.L0inv.T, 2)
    def step(self):
        g = self.rng.standard_normal(self.d)
        rem = 1.0
        while rem > 1e-9 and self.Qc.shape[1] < self.d:
            pg = g - self.Qc @ (self.Qc.T @ g)
            D = self.unvec(math.sqrt(self.dt) * rem * pg)
            if self.rel(self.A + D - self.A0) > self.thr:
                self.reenum()
            rs = self.rel(D); frac = 1.0
            if self.rel(self.A + D - self.A0) > self.thr:
                frac = 0.9 * self.thr / rs
                D = D * frac
            dq = quad(D, self.X)
            neg = dq < 0
            s = 2.0
            if neg.any():
                sv = (self.q[neg] - 1.0) / (-dq[neg])
                s = max(sv.min(), 0.0)
            if s <= 1.0:
                self.A = self.A + s * D
                self.t += s * frac * rem * self.dt
                rem *= (1 - s * frac)
                qn = self.q + s * dq
                new = qn <= 1 + 1e-9
                self.add_contacts(self.X[new])
                self.X = self.X[~new]
                self.fix_drift()
                self.q = quad(self.A, self.X)
                self.events.append((self.t, self.V(), len(self.contacts), self.Qc.shape[1]))
            else:
                self.A = self.A + D
                self.q = self.q + dq
                self.t += frac * rem * self.dt
                rem *= (1 - frac)

def make_lattice(n, kind, rng):
    kap = kappa(n)
    if kind == 'zn':
        return np.eye(n) * kap ** (1.0 / n)
    p = 100003
    M = np.eye(n)
    M[n - 1, :n - 1] = rng.integers(0, p, size=n - 1)
    M[n - 1, n - 1] = p
    lll(M); lll(M)
    return M * (kap / p) ** (1.0 / n)

def run_one(n, kind, rng, maxsteps=400000):
    B = make_lattice(n, kind, rng)
    dtf = float(os.environ.get('DTF', '1.0'))   # optional time-step factor for discretization-sensitivity checks
    dt = 1.25e-3 / n * dtf
    sim = Sim(n, B, rng, dt)
    a0 = sim.a0
    TK = 16 * math.log(n) / n ** 2
    V0 = sim.V(); VT_raw = None; VT_sc = None; pairsT = None
    steps = 0
    while sim.Qc.shape[1] < sim.d and steps < maxsteps:
        sim.step(); steps += 1
        if VT_raw is None and sim.t >= TK: VT_raw = sim.V()
        if VT_sc is None and sim.t >= TK * a0 ** 2: VT_sc = sim.V(); pairsT = len(sim.contacts)
    frozen = sim.Qc.shape[1] == sim.d
    sim.fix_drift()
    Vf = sim.V()
    if VT_raw is None: VT_raw = Vf
    if VT_sc is None: VT_sc = Vf; pairsT = len(sim.contacts)
    # independent final verification: fresh LLL + enumeration from the original basis
    Lc = np.linalg.cholesky(sim.A)
    Bc = B.copy() / sim.lam1; M = Lc.T @ Bc; lll(M, extras=(Bc,))
    Z = enum_points(M.T @ M, 1 + 1e-6)
    qq = quad(sim.A, Z @ Bc.T)
    minq = float(qq.min()) if len(qq) else None
    nmin = int((np.abs(qq - 1) < 1e-6).sum())
    return dict(n=n, kind=kind, dt=dt, lam1=sim.lam1, a0=a0, V0=V0, Vfinal=Vf, frozen=bool(frozen),
                tfreeze=sim.t, tfreeze_scaled=sim.t / a0 ** 2, TK=TK, VT_raw=VT_raw, VT_scaled=VT_sc,
                pairsT=pairsT, npairs_tracked=len(sim.contacts), rank=sim.Qc.shape[1],
                final_minq=minq, final_npairs=nmin, steps=steps, nenum=sim.nenum, warn=sim.warn,
                events=[(round(e[0], 8), round(e[1], 6), e[2], e[3]) for e in sim.events])

if __name__ == '__main__':
    n = int(sys.argv[1]); kind = sys.argv[2]; maxruns = int(sys.argv[3]); cap = float(sys.argv[4])
    seed = int(sys.argv[5]); out = sys.argv[6]
    rng = np.random.default_rng(seed)
    t0 = time.time(); k = 0
    with open(out, 'a') as f:
        while k < maxruns and time.time() - t0 < cap:
            try:
                r = run_one(n, kind, rng)
            except Exception as e:   # record numerical failures instead of crashing (they are counted, not hidden)
                open(out + '.failures', 'a').write(json.dumps(dict(n=n, kind=kind, error=repr(e))) + '\n'); k += 1
                continue
            f.write(json.dumps(r) + '\n'); f.flush(); k += 1
    print(n, kind, 'runs', k, 'time', round(time.time() - t0, 1))
