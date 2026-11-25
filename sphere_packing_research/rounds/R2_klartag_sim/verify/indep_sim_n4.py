#!/usr/bin/env python3
"""Verifier's independent, from-scratch simulation of dA = pi_F(dW) on L = Z^n (lambda_1 = 1), A_0 = 1.001 Id.
Differences from R2's code: fixed candidate box |z_i| <= R (validated a posteriori by own Fincke-Pohst), projection via SVD,
and after each contact event the rest of the Gaussian increment is DISCARDED (fresh independent increment next), i.e. a
different Euler scheme with the same continuum limit.  Usage: indep_sim_n4.py n nruns seed dtfactor"""
import sys, math, json, itertools, time
import numpy as np
import myenum as E
n = int(sys.argv[1]); nruns = int(sys.argv[2]); seed = int(sys.argv[3]); dtf = float(sys.argv[4]); Rb = 3 if n >= 5 else 4
d = n * (n + 1) // 2; iu = np.triu_indices(n, 1)
Z = np.array([z for z in itertools.product(range(-Rb, Rb + 1), repeat=n) if any(z) and next(a for a in z if a != 0) > 0], dtype=float)
Vc = np.concatenate([Z * Z, math.sqrt(2) * Z[:, iu[0]] * Z[:, iu[1]]], axis=1)       # <A, z z^T> = Vc @ vec(A)
def mat(a):
    A = np.zeros((n, n)); A[iu] = a[n:] / math.sqrt(2); A = A + A.T; A[np.diag_indices(n)] = a[:n]; return A
rng = np.random.default_rng(seed); dt = 1.25e-3 / n * dtf; kap = E.kappa(n)
out = open(f'indep_n{n}_seed{seed}_dtf{dtf}.jsonl', 'w'); t0 = time.time()
for run in range(nruns):
    a = np.zeros(d); a[:n] = 1.001
    contact = np.zeros(len(Z), bool); Q = np.zeros((d, 0)); t = 0.0; q = Vc @ a; nev = 0
    while Q.shape[1] < d:
        g = rng.standard_normal(d)
        pg = g - Q @ (Q.T @ g)
        da = math.sqrt(dt) * pg
        dq = Vc @ da
        cand = (~contact) & (dq < 0)
        s = 2.0
        if cand.any():
            sv = (q[cand] - 1.0) / (-dq[cand]); s = max(float(sv.min()), 0.0)
        if s <= 1.0:
            a = a + s * da; t += s * dt; q = Vc @ a; nev += 1
            contact |= (q <= 1 + 1e-10)
            C = Vc[contact]
            U, sg, Vt = np.linalg.svd(C, full_matrices=False)
            r = int((sg > 1e-8 * sg[0]).sum()); Q = Vt[:r].T
            a = a + np.linalg.lstsq(C, 1.0 - C @ a, rcond=None)[0]; q = Vc @ a
        else:
            a = a + da; t += dt; q = q + dq
    A = mat(a)
    Zm, qm = E.fp_dfs(A, 1 + 1e-6)          # rigorous check over the whole lattice Z^n
    minq = float(qm.min()) if len(qm) else None
    V = kap / math.sqrt(np.linalg.det(A))
    out.write(json.dumps(dict(V=V, t=t, minq=minq, npairs=int((np.abs(qm - 1) < 1e-6).sum()), tracked=int(contact.sum()), nev=nev)) + '\n'); out.flush()
print('done', nruns, 'runs', round(time.time() - t0, 1), 's')
