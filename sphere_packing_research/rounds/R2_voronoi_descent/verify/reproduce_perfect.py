"""Re-run the R2 code with the original seeds for the first k samples, capture the PERFECT-stage form
(not saved in the raw data) and verify it with the independent enumerator; also checks reproducibility."""
import sys, os, json, math, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
import voronoi as V
from vlib import enum_bfs, kappa
cd = {2: 1 / (2 * 3 ** .5), 3: 1 / (4 * 2 ** .5), 4: 1 / 8, 5: 1 / (8 * 2 ** .5), 6: 1 / (8 * 3 ** .5), 7: 1 / 16, 8: 1 / 16,
      9: 1 / (16 * 2 ** .5), 10: 1 / (16 * 3 ** .5)}
Vb = {n: 2 ** n * kappa(n) * cd[n] for n in cd}
jobs = [(int(a.split(':')[0]), a.split(':')[1], int(a.split(':')[2]), int(a.split(':')[3])) for a in sys.argv[1:]]
res = []
for n, mode, seed, k in jobs:
    saved = json.load(open(os.path.join(here, '..', 'data', f'runs_n{n}_{mode}_s{seed}.json')))['recs']
    rng = np.random.default_rng(1000 * n + seed)
    N = n * (n + 1) // 2
    for i in range(min(k, len(saved))):
        A, S = V.random_start(n, rng)
        V0 = kappa(n) / math.sqrt(np.linalg.det(A))
        # independent check of the start: min must be 1
        X0, v0 = enum_bfs(A, 1 + 1e-7)
        A, S, steps = V.to_perfect(A, S, rng, mode, [])
        Ap = A.copy()
        Vp = kappa(n) / math.sqrt(np.linalg.det(Ap))
        X, vals = enum_bfs(Ap, 1 + 1e-7)
        pairs = len(X) // 2
        iu = np.triu_indices(n)
        rk = int(np.linalg.matrix_rank((X[:, iu[0]] * X[:, iu[1]]).astype(float)))
        lpv, _, _ = V.descent_lp(A, S)
        A, S, nst, nonedge, val, Vpath = V.greedy_descent(A, S, rng)
        Ve = kappa(n) / math.sqrt(np.linalg.det(A))
        V.eutaxy_lp(A, S)
        r = saved[i]
        rec = dict(n=n, mode=mode, seed=seed, i=i, start_min=float(v0.min()), V0_match=bool(abs(V0 - r['V_start']) < 1e-9 * V0),
                   Vp=Vp, Vp_saved=r['V_perfect'], Vp_match=bool(abs(Vp - r['V_perfect']) < 1e-8 * Vp),
                   perfect_min=float(vals.min()), perfect_pairs=pairs, pairs_saved=r['pairs_perfect'], rank=rk, N=N,
                   Ve=Ve, Ve_saved=r['V_ext'], Ve_match=bool(abs(Ve - r['V_ext']) < 1e-8 * Ve), Vp_le_best=bool(Vp <= Vb[n] * (1 + 1e-8)))
        res.append(rec)
        print(rec, flush=True)
tag = '_'.join(sys.argv[1:]).replace(':', '-')
json.dump(res, open(os.path.join(here, f'reproduce_{tag}.json'), 'w'))
print('SUMMARY', len(res), 'samples; reproduced V_perfect', sum(r['Vp_match'] for r in res), 'V_ext', sum(r['Ve_match'] for r in res),
      '; perfect-stage min in [%.10f, %.10f]' % (min(r['perfect_min'] for r in res), max(r['perfect_min'] for r in res)),
      '; pairs match', sum(r['perfect_pairs'] == r['pairs_saved'] for r in res), '; rank==N', sum(r['rank'] == r['N'] for r in res),
      '; Vp<=best', sum(r['Vp_le_best'] for r in res))
