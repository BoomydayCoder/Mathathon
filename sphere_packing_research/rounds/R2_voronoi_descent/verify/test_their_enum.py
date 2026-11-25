"""Direct adversarial test of the R2 enumerator (LLL + C Fincke-Pohst) against the independent BFS enumerator
on random forms (well- and ill-conditioned, incl. forms close to the boundary of the PD cone, non-reduced bases)."""
import sys, os, numpy as np, time
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
import voronoi as V
from vlib import enum_bfs
rng = np.random.default_rng(2024)
bad = 0; tot = 0; maxcnt = 0; t0 = time.time()
for n in range(2, 11):
    for trial in range(60):
        kind = trial % 4
        G = rng.standard_normal((n, n))
        A = G.T @ G
        if kind == 1:   # ill conditioned
            A = A + 1e-3 * np.eye(n)
            D = np.diag(np.exp(rng.uniform(-2, 2, n))); A = D @ A @ D
        if kind == 2:   # random unimodular transformation of a nice form (non reduced basis)
            U = np.eye(n, dtype=np.int64)
            for _ in range(3 * n):
                i, j = rng.choice(n, 2, replace=False); U[:, i] += int(rng.integers(-2, 3)) * U[:, j]
            A0 = np.eye(n) + 0.5 * (np.ones((n, n)) - np.eye(n)) * (rng.random() < 0.5)
            A = U.T @ A0 @ U
        if kind == 3:   # nearly singular: A + t B close to tmax
            B = rng.standard_normal((n, n)); B = (B + B.T) / 2
            Li = np.linalg.inv(np.linalg.cholesky(A)); w = np.linalg.eigvalsh(Li @ B @ Li.T)
            if w[0] < 0: A = A + (-1 / w[0]) * 0.97 * B
        A = (A + A.T) / 2
        # scale so that the bound is a bit above the minimum
        Ar, U = V.lll_gram(A)
        b = np.diag(Ar).min()
        for fac in (1.0 + 1e-9, 1.3, 2.0):
            bound = b * fac
            X1, v1, fl = V.short_vectors(A, bound, maxc=200000)
            X2, v2 = enum_bfs(Ar, bound, tol=1e-9)
            tot += 1
            # compare with tolerance at the boundary: count vectors clearly inside
            c1 = int((v1 <= bound * (1 - 1e-7)).sum()); c2 = int((v2 <= bound * (1 - 1e-7)).sum())
            m1 = v1.min() if len(v1) else None; m2 = v2.min() if len(v2) else None
            maxcnt = max(maxcnt, c1)
            if 2 * c1 != c2 or (m1 is not None and abs(m1 - m2) > 1e-8 * b):
                bad += 1
                print('MISMATCH n', n, 'kind', kind, 'fac', fac, 'theirs', c1, 'mine/2', c2 / 2, m1, m2)
print('tests', tot, 'mismatches', bad, 'largest count', maxcnt, 'time', round(time.time() - t0, 1))
