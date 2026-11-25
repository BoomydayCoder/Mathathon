"""Fresh-seed replication of the R2 pipeline (their code, new seeds) for n=9,10; every perfect-stage and final form is
re-checked with the independent enumerator."""
import sys, os, json, math, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
import voronoi as V
from vlib import kappa, enum_bfs
out = []
for n, mode, budget in ((9, 'descent', 36), (9, 'random', 36), (10, 'descent', 40), (10, 'random', 40)):
    rng = np.random.default_rng(424200 + 10 * n + (mode == 'random'))
    t0 = time.time()
    while time.time() - t0 < budget:
        A, S = V.random_start(n, rng)
        A, S, steps = V.to_perfect(A, S, rng, mode, None)
        X, vals = enum_bfs(A, 1 + 1e-7)
        Vp = kappa(n) / math.sqrt(np.linalg.det(A)); okp = abs(vals.min() - 1) < 1e-7 and len(X) == 2 * len(S)
        A, S, nst, nonedge, val, Vpath = V.greedy_descent(A, S, rng)
        X, vals = enum_bfs(A, 1 + 1e-7)
        Ve = kappa(n) / math.sqrt(np.linalg.det(A)); oke = abs(vals.min() - 1) < 1e-7 and len(X) == 2 * len(S)
        out.append(dict(n=n, mode=mode, Vp=Vp, Ve=Ve, pairs_p=None, pairs_e=len(S), steps=nst, okp=bool(okp), oke=bool(oke)))
    R = [o for o in out if o['n'] == n and o['mode'] == mode]
    print(n, mode, 'runs', len(R), 'Vp median %.3f mean %.3f max %.3f' % (np.median([o['Vp'] for o in R]), np.mean([o['Vp'] for o in R]), max(o['Vp'] for o in R)),
          '| Ve median %.3f mean %.3f max %.3f' % (np.median([o['Ve'] for o in R]), np.mean([o['Ve'] for o in R]), max(o['Ve'] for o in R)),
          '| mean LP steps %.2f' % np.mean([o['steps'] for o in R]), '| independent min/pairs check ok:', all(o['okp'] and o['oke'] for o in R), flush=True)
    json.dump(out, open(os.path.join(here, 'fresh_control.json'), 'w'))
