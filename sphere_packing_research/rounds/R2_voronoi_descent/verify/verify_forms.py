"""Independent re-verification of ALL saved final forms A_ext (own enumeration, exact integer arithmetic)."""
import numpy as np, json, glob, math, os, sys, time
from scipy.optimize import linprog
from vlib import *
here = os.path.dirname(os.path.abspath(__file__))
delta = {2: 1 / (2 * math.sqrt(3)), 3: 1 / (4 * math.sqrt(2)), 4: 1 / 8, 5: 1 / (8 * math.sqrt(2)),
         6: 1 / (8 * math.sqrt(3)), 7: 1 / 16, 8: 1 / 16, 9: 1 / (16 * math.sqrt(2)), 10: 1 / (16 * math.sqrt(3))}
Vbest = {n: 2 ** n * kappa(n) * d for n, d in delta.items()}
out = []
t0 = time.time()
for f in sorted(glob.glob(os.path.join(here, '..', 'data', 'runs_n*_*.json'))):
    d = json.load(open(f))
    n = d['n']; N = n * (n + 1) // 2
    for i, r in enumerate(d['recs']):
        A = np.array(r['A_ext'], dtype=float)
        res = dict(file=os.path.basename(f), idx=i, n=n, mode=d['mode'])
        res['sym_err'] = float(np.abs(A - A.T).max())
        q = find_denominator(A)
        res['q'] = q; res['q_stored'] = r['denom_ext']
        X, vals = enum_bfs(A, 1 + 1e-5)
        res['min_float'] = float(vals.min())
        res['n_le_bound'] = int(len(X))
        if q > 0:
            Z = np.rint(q * A).astype(np.int64)
            res['rat_err'] = float(np.abs(q * A - Z).max())
            vi = np.einsum('ij,jk,ik->i', X, Z, X)
            res['min_int'] = int(vi.min())
            S = X[vi == vi.min()]
            dZ = int_det(Z.tolist())
            res['detZ'] = dZ
            detA = dZ / q ** n
            # normalise by the exact minimum (should be q)
            res['V_exact'] = kappa(n) / math.sqrt(dZ / float(vi.min()) ** n)
        else:
            S = X[vals <= vals.min() + 1e-7]
            res['V_exact'] = kappa(n) / math.sqrt(np.linalg.det(A / vals.min()))
        res['V_float'] = kappa(n) / math.sqrt(np.linalg.det(A / vals.min()))
        # one of each +- pair
        key = [tuple(x) for x in S]
        half = [x for x in key if x > tuple(-np.array(x))]
        Sh = np.array(half, dtype=np.int64)
        res['pairs'] = len(half); res['pairs_stored'] = r['pairs_ext']
        res['V_stored'] = r['V_ext']
        res['V_relerr'] = abs(res['V_exact'] - r['V_ext']) / r['V_ext']
        iu = np.triu_indices(n)
        Phi = (Sh[:, iu[0]] * Sh[:, iu[1]]).astype(float)
        res['rank'] = int(np.linalg.matrix_rank(Phi)); res['N'] = N
        # eutaxy LP: max t, A^-1 = sum lam_x x x^T, lam_x >= t
        m = len(Sh)
        Ainv = np.linalg.inv(A)
        Aeq = np.hstack([Phi.T, np.zeros((N, 1))])
        c = np.zeros(m + 1); c[-1] = -1
        Aub = np.hstack([-np.eye(m), np.ones((m, 1))])
        lp = linprog(c, A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=Ainv[iu], bounds=[(None, None)] * (m + 1), method='highs')
        res['eut'] = float(lp.x[-1]) if lp.status == 0 else float('nan')
        res['eut_stored'] = r['eut_lp']
        res['index'] = int(lattice_index(Sh, n))
        res['maxcoord'] = int(np.abs(Sh).max())
        if q > 0:
            ip = Sh @ Z @ Sh.T
            res['ipset'] = sorted(set(int(a) for a in np.unique(np.abs(ip))))
        res['exceeds_best'] = bool(res['V_exact'] > Vbest[n] * (1 + 1e-9))
        res['V_over_best'] = res['V_exact'] / Vbest[n]
        res['V_perfect'] = r['V_perfect']; res['descent_steps'] = r['descent_steps']
        out.append(res)
    print(os.path.basename(f), 'done', round(time.time() - t0, 1), flush=True)
json.dump(out, open(os.path.join(here, 'verify_forms_results.json'), 'w'))
# ---- report
import collections
print('\nTOTAL records verified:', len(out))
for n in range(2, 11):
    R = [o for o in out if o['n'] == n]
    if not R: continue
    print(f"n={n}: {len(R)} forms | min of (exact min / q) deviation: min_int==q in {sum(o.get('min_int')==o['q'] for o in R)}/{len(R)}"
          f" | min_float range [{min(o['min_float'] for o in R):.10f},{max(o['min_float'] for o in R):.10f}]"
          f" | max rat_err {max(o.get('rat_err',0) for o in R):.2e}"
          f" | pairs match {sum(o['pairs']==o['pairs_stored'] for o in R)}/{len(R)}"
          f" | q match {sum(o['q']==o['q_stored'] for o in R)}/{len(R)}"
          f" | max V relerr {max(o['V_relerr'] for o in R):.2e}"
          f" | perfect(rank=N) {sum(o['rank']==o['N'] for o in R)}/{len(R)}"
          f" | eutactic(eut>1e-9) {sum(o['eut']>1e-9 for o in R)}/{len(R)} min eut {min(o['eut'] for o in R):.3e}"
          f" | index==1 {sum(o['index']==1 for o in R)}/{len(R)} | maxcoord {max(o['maxcoord'] for o in R)}"
          f" | max V {max(o['V_exact'] for o in R):.5f} (best {Vbest[n]:.5f}) exceed: {sum(o['exceeds_best'] for o in R)}"
          f" | q values {sorted(collections.Counter(o['q'] for o in R).items())}")
