#!/usr/bin/env python3
# Referee: independent exhaustive t-grid enumeration on a sublattice (rows of F in A's saved basis),
# using the R2 LLL + Fincke-Pohst code (not the builders' enumeration code).
import sys, json, time, numpy as np
from math import gcd
sys.path.insert(0, '/home/user/Mathathon/sphere_packing_research/rounds/R2_voronoi_descent')
import voronoi
R = '/home/user/Mathathon/sphere_packing_research/rounds/R10_weil_lattices/'
p = int(sys.argv[1]); Ff = sys.argv[2]; step = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
tl = float(sys.argv[4]) if len(sys.argv) > 4 else 400
t0 = time.time()
B = np.load(R + 'A/B_p%d_L0.npy' % p); E = np.load(R + 'A/E_p%d_L0.npy' % p); F = np.load(R + 'A/' + Ff)
Dtab = {5: 5, 13: 117, 17: 153, 29: 1421, 37: 2997, 61: 13725}; D = Dtab[p]   # verified exactly in ref.py
eps = {5: (1 + 5 ** .5) / 2, 13: (3 + 13 ** .5) / 2, 17: 4 + 17 ** .5, 29: (5 + 29 ** .5) / 2, 37: 6 + 37 ** .5, 61: (39 + 5 * 61 ** .5) / 2}[p]
THR = (p - 1) ** 2 * D; theta = (p - 1) ** 2 / (4 * p)
Bo = B.astype(object); Eo = E.astype(object)
k = F.shape[0]
# primitivity of F (saturation): gcd of several k x k minors
def bareiss(M):
    nn = len(M); A = [list(map(int, r)) for r in M]; prev = 1; s = 1
    for q in range(nn):
        if A[q][q] == 0:
            for i in range(q + 1, nn):
                if A[i][q] != 0: A[q], A[i] = A[i], A[q]; s = -s; break
            else: return 0
        akk = A[q][q]; Ak = A[q]
        for i in range(q + 1, nn):
            Ai = A[i]; aik = Ai[q]
            for j in range(q + 1, nn): Ai[j] = (Ai[j] * akk - aik * Ak[j]) // prev
        prev = akk
    return s * A[nn - 1][nn - 1]
rng = np.random.default_rng(7); g = 0
for trial in range(6):
    cols = rng.choice(F.shape[1], size=k, replace=False); g = gcd(g, abs(bareiss(F[:, cols].tolist())))
    if g == 1: break
print('F', Ff, F.shape, 'gcd of sampled maximal minors =', g, '(1 => primitive/saturated sublattice)', flush=True)
Fo = F.astype(object)
GB = np.array(Fo.dot(Bo).dot(Fo.T), dtype=float); GE = np.array(Fo.dot(Eo).dot(Fo.T), dtype=float)
Gp = (GB + GE / np.sqrt(p)) / 2; Gm = (GB - GE / np.sqrt(p)) / 2
bound = 2 * np.sqrt(theta * D) * np.cosh(step / 2) * (1 + 1e-9)
L = np.log(eps); ts = np.arange(-L, L + step / 2 + 1e-12, step)
print('p', p, 'period 2 log eps =', 2 * L, 'grid', len(ts), 'points from', ts[0], 'to', ts[-1], 'bound', bound, flush=True)
seen = {}; minN = None; tot = 0; done = 0
for t in ts:
    if time.time() - t0 > tl: print('TIME LIMIT reached at t =', t); break
    G = np.exp(t) * Gp + np.exp(-t) * Gm
    try:
        X, v, flag = voronoi.short_vectors(G, bound, maxc=20000, maxnodes=2_000_000_000)
    except RuntimeError as ex:
        print('t=%.3f node limit' % t); continue
    if flag: print('t=%.3f flag' % t, flag)
    done += 1; tot += len(X)
    for x in X:
        w = (x.astype(object)).dot(Fo); b = int(w.dot(Bo).dot(w)); e = int(w.dot(Eo).dot(w)); N = p * b * b - e * e
        key = tuple(int(c) for c in w)
        if key[next(i for i, c in enumerate(key) if c != 0)] < 0: key = tuple(-c for c in key)
        seen[key] = (N, b, e)
        minN = N if minN is None else min(minN, N)
viol = {kk: vv for kk, vv in seen.items() if vv[0] < THR}
from collections import Counter
print('SUMMARY p=%d %s rank %d: t-points completed %d/%d, vectors (with mult.) %d, distinct +-pairs %d, min N = %s, THR = %d, min Gamma/theta = %s, distinct violating +-pairs = %d, (N,b,e) histogram: %s, time %.1fs, nodes %d'
      % (p, Ff, k, done, len(ts), tot, len(seen), minN, THR, (minN / THR if minN else None), len(viol), dict(Counter(seen.values()).most_common(8)), time.time() - t0, voronoi.STATS['nodes']), flush=True)
