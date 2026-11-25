"""p = 29 (r = 7, threshold 49; (p+1)/2 = 15 composite, known unstable): sample of short vectors of the full rank-56
lattice from the formula (NOT exhaustive: enumeration stopped after maxc vectors) with exact (norm, |s|, N(nu))."""
import sys, collections, time
sys.path.insert(0,'.'); sys.path.insert(0,'../R10_weil_lattices/orchestrator')
import ideal_model as IM, scan_exact as SE, numpy as np
p = 29; N = 28
Gi, Ei, Bz, F = IM.build(p, N, verbose=True)
U0 = SE.reduce_exact(Gi); Gr = U0.T @ Gi @ U0; Er = U0.T @ Ei @ U0
print('reduced diag:', sorted(collections.Counter(int(Gr[i, i]) for i in range(2 * N)).items()), flush=True)
for bound in (6.01,):
    t0 = time.time()
    try:
        Xs, vals, flag = IM.V.short_vectors(np.array(Gr, dtype=float), bound, maxc=100000, maxnodes=30_000_000_000)
    except Exception as ex:
        print('bound', bound, 'failed', ex, flush=True); continue
    cnt = collections.Counter()
    for x in np.array(Xs, dtype=object):
        nrm = int(x @ Gr @ x); s = int(x @ Er @ x) // 2
        cnt[(nrm, abs(s), (p * nrm * nrm - 4 * s * s) // 16)] += 1
    print('bound %.2f: %d vectors (flag %d: 2 = truncated sample), (norm,|s|,N(nu)): %s  [%.0fs]' % (bound, len(Xs), flag, dict(sorted(cnt.items())), time.time() - t0), flush=True)
