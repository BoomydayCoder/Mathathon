#!/usr/bin/env python3
"""LP-only scan (trace-LP vertex of the Ryshkov polyhedron for random lattices), independent seeds, more dimensions.
Usage: python3 lp_scan.py n kind nsamp timecap seed outfile"""
import sys, os, json, time, math
import numpy as np
import tracelp as T
from tracelp import kappa, trace_lp, verify, describe, gram_random, F_basis
# Conway-Sloane table of best known lattice centre densities (laminated lattices) for n = 15..20
T.DELTA.update({15: 1 / (16 * math.sqrt(2)), 16: 1 / 16, 17: 1 / 16, 18: 1 / (8 * math.sqrt(3)), 19: 1 / (8 * math.sqrt(2)), 20: 1 / 8,
                21: 1 / (4 * math.sqrt(2)), 22: 1 / (2 * math.sqrt(3)), 23: 0.5, 24: 1.0,
                26: 1 / math.sqrt(3), 28: 2.0 / 3})   # Lambda_26, Lambda_28 (conservative: true lattice densities)
n = int(sys.argv[1]); kind = sys.argv[2]; nsamp = int(sys.argv[3]); cap = float(sys.argv[4]); seed = int(sys.argv[5]); out = sys.argv[6]
rng = np.random.default_rng(100000 * n + 7 * seed + (1 if kind == 'gm' else 2) + 5000)
N = n * (n + 1) // 2; Vb = T.Vbest(n)
GMP = int(os.environ.get('GMP', '100003'))     # modulus of the Goldstein-Mayer family (sensitivity check)
label = kind if GMP == 100003 else f'{kind}_p{GMP}'
t0 = time.time(); recs = []
for i in range(nsamp):
    if time.time() - t0 > cap:
        break
    Q = gram_random(n, kind, rng, p=GMP)
    ts = time.time()
    rec = dict(n=n, kind=kind, idx=i, Q=np.round(Q, 12).tolist())
    try:
        A, S, lam2, info = trace_lp(Q)
        mn, npairs = verify(A)
        d = describe(A, Q, n)
        rec.update(lam1sq=lam2, V0=kappa(n) * lam2 ** (n / 2), V_LP=d['V'], v_LP=d['v'], B0sq=d['B0sq'],
                   eig_ratio=d['eig_ratio'], pairs_LP=npairs, rank_LP=int(F_basis(S)[1]), min_LP=mn,
                   lp_info=info, A_LP=np.round(A, 12).tolist(), t_lp=time.time() - ts)
        assert abs(mn - 1) < 1e-9 and d['V'] <= Vb * (1 + 1e-7), ('verification failed', mn, d['V'], Vb)
    except Exception as e:
        rec['lp_error'] = repr(e)
        print('LP FAIL', n, kind, i, repr(e), flush=True)
    recs.append(rec)
    json.dump(dict(n=n, kind=label, seed=seed, N=N, Vbest=Vb, elapsed=time.time() - t0, recs=recs), open(out, 'w'))
ok = [r for r in recs if 'V_LP' in r]
print(f'lpscan n={n} {label}: {len(ok)}/{len(recs)} LPs in {time.time()-t0:.1f}s  V_LP/N mean={np.mean([r["V_LP"] for r in ok])/N:.3f} '
      f'pairs/N={np.mean([r["pairs_LP"] for r in ok])/N:.3f} Vbest/N={Vb/N:.3f}', flush=True)
