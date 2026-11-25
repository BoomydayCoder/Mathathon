"""Validation of trace_lp on known cases: if Q is an extreme (perfect + eutactic) form, Q^-1 lies in the interior of
its own Voronoi domain, so the trace-LP optimum must be Q/min(Q) itself."""
import numpy as np, math, json, sys
import tracelp as T
from tracelp import kappa, trace_lp, verify, describe, Vbest

def An(n):
    return 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
def Dn(n):
    G = An(n); G[n - 1, n - 2] = G[n - 2, n - 1] = 0; G[n - 1, n - 3] = G[n - 3, n - 1] = -1; return G
def En(n):
    # E_n Dynkin: chain 1..n-1, node n attached to node 3
    G = 2 * np.eye(n)
    for i in range(n - 2):
        G[i, i + 1] = G[i + 1, i] = -1
    G[n - 1, 2] = G[2, n - 1] = -1
    return G
out = []
for name, G, Vexp, pexp in [('A6', An(6), kappa(6) * 2 ** 3 / math.sqrt(7), 21), ('D6', Dn(6), kappa(6) * 2 ** 3 / 2, 30),
                            ('E6', En(6), Vbest(6), 36), ('E7', En(7), Vbest(7), 63), ('E8', En(8), Vbest(8), 120),
                            ('D4', Dn(4), kappa(4) * 4 / 2, 12), ('A10', An(10), kappa(10) * 2 ** 5 / math.sqrt(11), 55)]:
    n = G.shape[0]
    Q = G / np.linalg.det(G) ** (1.0 / n)
    A, S, lam2, info = trace_lp(Q)
    mn, npairs = verify(A)
    d = describe(A, Q, n)
    ok = abs(d['V'] - Vexp) < 1e-6 * Vexp and npairs == pexp and abs(mn - 1) < 1e-9 and d['B0sq'] < 1e-12
    out.append(dict(name=name, V=d['V'], Vexp=Vexp, pairs=npairs, pairs_exp=pexp, minA=mn, B0sq=d['B0sq'], rounds=info['rounds'], ok=bool(ok)))
    print(name, 'V=%.6f expected %.6f pairs=%d (exp %d) min=%.12f B0sq=%.2e rounds=%d ncons=%d  %s' % (
        d['V'], Vexp, npairs, pexp, mn, d['B0sq'], info['rounds'], info['ncons'], 'OK' if ok else 'MISMATCH'))
json.dump(out, open('data/validation.json', 'w'), indent=1)
