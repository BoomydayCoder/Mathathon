"""p = 29 (threshold r^2 = 49, known unstable; Round 10 reports vector violators with Gamma/theta = 5/7, i.e. N(nu) = 35):
targeted search in the full rank-56 lattice from the formula at the point of the t-circle where a (norm 6, s = +-11)
vector is shortest (t = -+0.83), bound 4.45.  Confirms that (A') fails at p = 29 in the model."""
import sys, collections, time, math
sys.path.insert(0,'.'); sys.path.insert(0,'../R10_weil_lattices/orchestrator')
import ideal_model as IM, scan_exact as SE, numpy as np
p = 29; N = 28
Gi, Ei, Bz, F = IM.build(p, N, verbose=True)
U0 = SE.reduce_exact(Gi)
for t in (-0.83, 0.83):
    t0 = time.time(); U = U0.copy()
    Gt = SE.qt_float(U.T @ Gi @ U, U.T @ Ei @ U, t, p)
    _, Uf = IM.V.lll_gram(Gt); U = U @ np.array(Uf, dtype=object)
    Gt = SE.qt_float(U.T @ Gi @ U, U.T @ Ei @ U, t, p)
    try:
        Xs, vals, flag = IM.V.short_vectors(Gt, 4.45, maxc=300000, maxnodes=40_000_000_000)
    except Exception as ex:
        print('t=%.2f failed: %r' % (t, ex), flush=True); continue
    cnt = collections.Counter()
    for x in np.array(Xs, dtype=object):
        w = U @ x; nrm = int(w @ Gi @ w); s = int(w @ Ei @ w) // 2
        cnt[(nrm, s, (p * nrm * nrm - 4 * s * s) // 16)] += 1
    print('t=%+.2f: %d vectors with Q_t <= 4.45 (flag %d); (norm, s, N(nu)): %s  [%.0fs]' % (t, len(Xs), flag, dict(sorted(cnt.items())), time.time() - t0), flush=True)
