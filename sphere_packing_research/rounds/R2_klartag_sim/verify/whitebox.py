#!/usr/bin/env python3
"""White-box tests of R2's sim.py building blocks: Brownian covariance on Sym_n, isometry vec<->Tr, projection, enumeration."""
import sys, math, json
import numpy as np
sys.dont_write_bytecode = True
import sim_r2_copy as S
import myenum as E
rng = np.random.default_rng(12345)
n = 5; B = S.make_lattice(n, 'gm', rng); sm = S.Sim(n, B, rng, 1e-4)
d = sm.d
# 1. unvec/vecA are inverse and vecA is an isometry (Sym_n, Tr(AB)) -> R^d
M1 = rng.standard_normal((n, n)); M1 = M1 + M1.T; M2 = rng.standard_normal((n, n)); M2 = M2 + M2.T
print('isometry err', abs(sm.vecA(M1) @ sm.vecA(M2) - np.trace(M1 @ M2)), 'roundtrip err', np.abs(sm.unvec(sm.vecA(M1)) - M1).max())
# 2. covariance of unvec(g): diagonal variance 1, off-diagonal variance 1/2, E Tr(W^2) = d
G = rng.standard_normal((200000, d)); Ws = np.array([sm.unvec(g) for g in G[:50000]])
print('var diag', Ws[:, 0, 0].var(), Ws[:, 2, 2].var(), 'var offdiag', Ws[:, 0, 1].var(), Ws[:, 1, 3].var(), 'sym', np.abs(Ws - Ws.transpose(0, 2, 1)).max(),
      'E Tr W^2', np.mean(np.einsum('kij,kij->k', Ws, Ws)), 'd', d, 'cov(W00,W11)', np.mean(Ws[:, 0, 0] * Ws[:, 1, 1]), 'cov(W01,W02)', np.mean(Ws[:, 0, 1] * Ws[:, 0, 2]))
# 3. vec(x).vecA(S) = x^T S x
X = rng.standard_normal((7, n)); print('quad err', np.abs(sm.vec(X) @ sm.vecA(M1) - S.quad(M1, X)).max())
# 4. projection: after adding contacts, the projected increment annihilates all contact points and is the orthogonal projection
Xc = rng.integers(-3, 4, size=(6, n)).astype(float); sm.contacts = np.zeros((0, n)); sm.Qc = np.zeros((d, 0)); sm.add_contacts(Xc)
g = rng.standard_normal(d); pg = g - sm.Qc @ (sm.Qc.T @ g); Dm = sm.unvec(pg)
print('rank', sm.Qc.shape[1], 'max |D x.x| on contacts', np.abs(S.quad(Dm, Xc)).max(), 'orthonormal err', np.abs(sm.Qc.T @ sm.Qc - np.eye(sm.Qc.shape[1])).max(),
      'idempotent/orth: <pg, g-pg> =', pg @ (g - pg))
# explicit orthogonal projector from the definition F = {B: <B, x x^T> = 0}
Mc = np.array([sm.vecA(np.outer(x, x)) for x in Xc]); P = np.eye(d) - Mc.T @ np.linalg.pinv(Mc @ Mc.T) @ Mc
print('projection vs definition err', np.abs(P @ g - pg).max())
# 5. their enum_points versus my DFS and box enumeration on random forms/lattices (counts and min), incl. bound 2 (candidate lists)
bad = 0; tests = 0
for nn in (2, 3, 4, 5, 6, 7, 8):
    for trial in range(6):
        Bq = S.make_lattice(nn, 'gm', rng); T = rng.standard_normal((nn, nn)); A = T @ T.T / nn + 0.3 * np.eye(nn)
        Lc = np.linalg.cholesky(A); M = Lc.T @ Bq; Bc = Bq.copy(); S.lll(M, extras=(Bc,))
        Gm = M.T @ M
        lam = np.sort(S.quad(Gm, S.enum_points(Gm, Gm.diagonal().min())))[0]
        for bound in (lam * 1.0000001, lam * 2.0, lam * 3.5):
            Z1 = S.enum_points(Gm, bound); Z2, q2 = E.fp_dfs(Gm, bound)
            (bx, b, size) = E.box_enum(Gm, bound)
            s1 = set(map(tuple, np.round(Z1).astype(int) * np.sign(np.round(Z1)[np.arange(len(Z1)), [np.nonzero(r)[0][0] for r in np.round(Z1)]])[:, None].astype(int)))
            s2 = set(map(tuple, np.round(Z2).astype(int)))
            ok = (s1 == s2)
            if bx is not None:
                s3 = set(map(tuple, np.round(bx[0]).astype(int))); ok = ok and (s3 == s2)
            tests += 1; bad += (not ok)
print('enumeration cross-check: tests', tests, 'mismatches', bad)
