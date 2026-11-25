#!/usr/bin/env python3
# Task (4) heuristic on the FULL rank-120 even unimodular member (t=+log eps): G-orbits of the norm-10 vectors found in Fix(C_5),
# inner-product spectrum between orbit elements; any |ip| >= 6 gives a vector of norm <= 8 (pair difference).
import sys, os, json
from iso import *
p = 61; results = []
Wl = setup(p); rep = Rep(Wl); n = Wl.n; g0 = Wl.g; sp = math.sqrt(p); le = math.log(Wl.eps); sD = math.sqrt(Wl.D)
Bf = Wl.B.astype(float); Ef = Wl.E.astype(float)
Qf = (math.exp(le) * (Bf + Ef / sp) / 2 + math.exp(-le) * (Bf - Ef / sp) / 2) / sD
Qi = np.rint(Qf).astype(np.int64); assert np.abs(Qf - Qi).max() < 1e-6
def qs(Wv):
    b = np.einsum('ij,jk,ik->i', Wv, Wl.B, Wv).astype(float); e = np.einsum('ij,jk,ik->i', Wv, Wl.E, Wv).astype(float)
    mp_ = math.exp(le) * (b + e / sp) / 2 / sD; mm_ = math.exp(-le) * (b - e / sp) / 2 / sD
    return mp_ + mm_, (mp_ - mm_) * sp / 2, mp_ * mm_ / Wl.theta
split = cyclic_components(rep, (g0, 0, 0, pow(g0, p - 2, p)), 'split torus C_60', results, sums=True)
Et, N, r = split[(4, 12)]
C, exact = sat_basis(Wl, Et, N, r['dimU'], int(r['glue_index']))
Q = C.astype(float) @ Qf @ C.T.astype(float)
U, _ = lll_int(np.rint(Q * 2.0 ** 30).astype(np.int64)); Uf = U.astype(float)
X, vals, flag, nodes = fp(Uf @ Q @ Uf.T, 10 + 1e-6, maxc=400000)
V = X @ U @ C
q, s, gt = qs(V); log('norm-10 vectors in Fix(C_5): %d pairs; |s| values %s' % (len(V), sorted(set(np.rint(np.abs(s)).astype(int).tolist()))))
PSf = rep.PS.astype(float); PTf = rep.PT.astype(float)
def canon(v):
    nz = v[np.nonzero(v)[0][0]]
    return v if nz > 0 else -v
def orbit(x, cap=130000):
    x = canon(x); seen = {x.tobytes()}; allv = [x[None, :]]; fr = x[None, :].astype(float)
    while len(fr):
        Y = np.rint(np.vstack([fr @ PSf, fr @ PTf])).astype(np.int64)
        new = []
        for y in Y:
            y = canon(y); k = y.tobytes()
            if k not in seen: seen.add(k); new.append(y)
        if not new: break
        new = np.array(new); allv.append(new); fr = new.astype(float)
        assert len(seen) <= cap
    return np.vstack(allv), seen
keys_done = set(); orbs = []
T0_ = time.time()
for i in range(len(V)):
    k = canon(V[i]).tobytes()
    if k in keys_done: continue
    O, seen = orbit(V[i]); keys_done |= seen; orbs.append(O)
    qo, so, go = qs(O[:1])
    log('orbit %d: %d vectors (mod +-), rep norm %.3f |s| %.3f Gamma/theta %.4f' % (len(orbs), len(O), qo[0], abs(so[0]), go[0]))
    if len(orbs) >= 40 or time.time() - T0_ > 90: break
out = dict(orbits=[len(O) for O in orbs], pairs=[])
best = None
for a in range(len(orbs)):
    x = orbs[a][0]; qx = Qi @ x
    for b_ in range(a, len(orbs)):
        ip = np.abs(orbs[b_] @ qx)
        hist = {int(v): int(c) for v, c in zip(*np.unique(ip, return_counts=True))}
        log('ip spectrum orbit %d rep vs orbit %d: %s' % (a + 1, b_ + 1, hist))
        out['pairs'].append(dict(a=a + 1, b=b_ + 1, hist=hist))
        for val in sorted(hist, reverse=True):
            if val >= 10 or val < 5: continue
            idx = np.nonzero(ip == val)[0][:2000]
            sg = np.sign(orbs[b_][idx] @ qx)[:, None]
            Wd = x[None, :] - sg * orbs[b_][idx]
            qd, sd, gd = qs(Wd)
            j = int(np.argmin(gd))
            log('   |ip|=%d: differences have norm %.3f, |s| values %s, min Gamma/theta %.4f' % (val, qd[j], sorted(set(np.rint(np.abs(sd)).astype(int).tolist()))[:12], gd[j]))
            if best is None or gd[j] < best[0]: best = (float(gd[j]), float(qd[j]), float(abs(sd[j])), Wd[j].tolist())
out['best'] = best
json.dump(out, open(os.path.join(HERE, 'unimod_full_orbit_p61.json'), 'w'), indent=1)
log('BEST over pair differences: Gamma/theta %.4f norm %.2f |s| %.2f' % (best[0], best[1], best[2]) if best else 'no pair with 5 <= |ip| < 10')
log('DONE extra4')
