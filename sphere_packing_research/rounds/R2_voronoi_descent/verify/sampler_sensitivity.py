"""How sampler-specific are the 'lands on the best lattice' frequencies?
Variant 'R2'   : the R2 walk as coded (direction isotropic w.r.t. Tr in the CURRENT LLL-reduced Z^n coordinates), fresh seeds (control).
Variant 'fixed': Klartag-faithful ray walk: the lattice L = M Z^n sits in a FIXED Euclidean R^n (start form = Id on L, i.e. a_0 Id
                 touching L), directions isotropic w.r.t. Tr(B'^2) in those fixed Euclidean coordinates, projected on
                 F' = {B' : B'[y]=0 for contact points y in L}.  (Only the direction law differs; the line search is the R2 code.)
"""
import sys, os, json, math, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
import voronoi as V
from vlib import kappa, enum_bfs
cd = {2: 1 / (2 * 3 ** .5), 3: 1 / (4 * 2 ** .5), 4: 1 / 8, 5: 1 / (8 * 2 ** .5), 6: 1 / (8 * 3 ** .5), 7: 1 / 16, 8: 1 / 16,
      9: 1 / (16 * 2 ** .5), 10: 1 / (16 * 3 ** .5)}
Vb = {n: 2 ** n * kappa(n) * cd[n] for n in cd}
state = {}
_orig_rebase = V.rebase
def tracking_rebase(A, S):
    G, U = V.lll_gram(A)
    Ui = np.rint(np.linalg.inv(U)).astype(np.int64)
    if 'M' in state:
        state['M'] = state['M'] @ U
    return G, S @ Ui.T
V.rebase = tracking_rebase

def to_perfect_fixed(A, S, rng, mode):
    n = A.shape[0]; N = n * (n + 1) // 2
    state['M'] = np.linalg.cholesky(A).T.copy()      # M^T M = A : lattice M Z^n with Euclidean form Id
    while True:
        M = state['M']
        # consistency: Gram of current basis w.r.t. current Euclidean form is not Id any more; we only need F' and the map B = M^T B' M
        Y = S @ M.T
        Fb, r = V.F_basis(Y)
        if r == N:
            break
        g = rng.standard_normal(Fb.shape[0])
        Bp = (g @ Fb).reshape(n, n); Bp = (Bp + Bp.T) / 2
        B = M.T @ Bp @ M
        B = (B + B.T) / 2
        B /= np.linalg.norm(B)
        if mode == 'descent':
            if np.trace(np.linalg.solve(A, B)) > 0: B = -B
        else:
            Li = np.linalg.inv(np.linalg.cholesky(A))
            if np.linalg.eigvalsh(Li @ B @ Li.T)[0] >= -1e-12: B = -B
        A, S, t = V.line_step(A, B)
    del state['M']
    return A, S

out = []
tl = float(sys.argv[1]) if len(sys.argv) > 1 else 150
t00 = time.time()
plan = [(4, 300), (5, 300), (6, 200), (7, 150), (8, 60)]
for n, ns in plan:
    for mode in ('descent', 'random'):
        for variant in ('R2', 'fixed'):
            rng = np.random.default_rng(777000 + 10 * n + (mode == 'random') + 100 * (variant == 'fixed'))
            Vs = []; t0 = time.time(); bad = 0
            budget = tl * {4: .04, 5: .06, 6: .08, 7: .12, 8: .2}[n] / 1.0
            for i in range(ns):
                if time.time() - t0 > budget: break
                A, S = V.random_start(n, rng)
                if variant == 'R2':
                    A, S, _ = V.to_perfect(A, S, rng, mode)
                else:
                    A, S = to_perfect_fixed(A, S, rng, mode)
                # independent check of the min and perfection
                X, vals = enum_bfs(A, 1 + 1e-7)
                if abs(vals.min() - 1) > 1e-7: bad += 1
                Vs.append(kappa(n) / math.sqrt(np.linalg.det(A)))
            Vs = np.array(Vs)
            fb = float(np.mean(np.abs(Vs - Vb[n]) < 1e-5 * Vb[n]))
            se = math.sqrt(max(fb * (1 - fb), 1e-12) / len(Vs))
            rec = dict(n=n, mode=mode, variant=variant, runs=len(Vs), freq_best=fb, se=se, median=float(np.median(Vs)), mean=float(Vs.mean()),
                       max=float(Vs.max()), bad_min=bad, distinct=len(set(np.round(Vs, 5))))
            out.append(rec)
            print(rec, flush=True)
json.dump(out, open(os.path.join(here, 'sampler_sensitivity.json'), 'w'), indent=1)
