#!/usr/bin/env python3
# Task (4), restricted version: short vectors of the EVEN UNIMODULAR member (t = +-log eps) inside the structured sublattices (rank <= 40).
# q = unimodular norm, s = (m~_+ - m~_-) sqrt(p)/2 ; Gamma = (q^2 - 4 s^2/p)/4 ; violator iff Gamma < theta.
import sys, os, json
from iso import *
p = int(sys.argv[1]); results = []; out = []
Wl = setup(p); rep = Rep(Wl); n = Wl.n; g0 = Wl.g; sp = math.sqrt(p); le = math.log(Wl.eps); sD = math.sqrt(Wl.D)
Bf = Wl.B.astype(float); Ef = Wl.E.astype(float)
Qfull = (math.exp(le) * (Bf + Ef / sp) / 2 + math.exp(-le) * (Bf - Ef / sp) / 2) / sD
log('unimodular member at t=+log eps: integrality residual %.2e, even diagonal: %s, float logdet %.2e' % (np.abs(Qfull - np.rint(Qfull)).max(), bool(np.all(np.rint(Qfull.diagonal()) % 2 == 0)), np.linalg.slogdet(Qfull)[1]))
log('   diag of unimodular Gram in saved (B-reduced) basis: min %.0f max %.0f' % (Qfull.diagonal().min(), Qfull.diagonal().max()))
todo = []
gl = find_gens(p, '2T'); o = group_components(rep, gl[0], 'SL_2(3)', CH_2T, results, expect_order=24)
todo += [(r['label'], Et, N, r) for (Et, N, r) in o.values() if r['dimU'] <= 40]
if p % 10 in (1, 9):
    gl = find_gens(p, '2I'); o = group_components(rep, gl[0], 'SL_2(5)', CH_2I, results, expect_order=120)
    todo += [(r['label'], Et, N, r) for (Et, N, r) in o.values() if r['dimU'] <= 40]
split = cyclic_components(rep, (g0, 0, 0, pow(g0, p - 2, p)), 'split torus C_%d' % (p - 1), results, sums=True)
for c, (Et, N, r) in split.items():
    if r and r['dimU'] <= 40 and c in ((4,), (4, 12), (4, 20), (12,)): todo.append((r['label'], Et, N, r))
todo.sort(key=lambda x: x[3]['dimU'])
for lab, Et, N, r in todo:
    k = r['dimU']
    C, exact = sat_basis(Wl, Et, N, k, int(r['glue_index']))
    Bs, Es = Wl.grams(C); Bs = Bs.astype(float); Es = Es.astype(float)
    for sgn in (1, -1):
        t = sgn * le
        Q = (math.exp(t) * (Bs + Es / sp) / 2 + math.exp(-t) * (Bs - Es / sp) / 2) / sD
        U, _ = lll_int(np.rint(Q * 2.0 ** 30).astype(np.int64)); Uf = U.astype(float); Q2 = Uf @ Q @ Uf.T
        found = None
        for bound in (8, 10, 12, 14):
            X, vals, flag, nodes = fp(Q2, bound + 1e-6, maxc=400000, maxnodes=3 * 10**8)
            if flag: log('   enumeration overflow at bound %d (%s)' % (bound, lab)); break
            if len(X):
                W_ = X @ U @ C
                b = np.einsum('ij,jk,ik->i', W_, Wl.B, W_).astype(float); e = np.einsum('ij,jk,ik->i', W_, Wl.E, W_).astype(float)
                mp_ = math.exp(t) * (b + e / sp) / 2 / sD; mm_ = math.exp(-t) * (b - e / sp) / 2 / sD
                q = mp_ + mm_; s = (mp_ - mm_) * sp / 2; gam = mp_ * mm_
                qm = q.min(); sel = np.abs(q - qm) < 1e-6
                hist = {}
                for qq, ss, gg in zip(np.rint(q).astype(int), np.rint(s * 2) / 2, gam):
                    key = (int(qq), abs(float(ss))); hist[key] = hist.get(key, 0) + 1
                found = dict(label=lab, rank=k, t_sign=sgn, bound=bound, min_norm=float(qm), n_pairs_at_bound=len(X), nodes=nodes,
                             hist_q_abs_s={str(kk): v for kk, v in sorted(hist.items())}, min_gamma_over_theta=float(gam.min() / Wl.theta),
                             s_integrality=float(np.abs(s - np.rint(s)).max()))
                break
        if found is None: found = dict(label=lab, rank=k, t_sign=sgn, note='no vector of unimodular norm <= %d found (or overflow)' % bound)
        log('UNIMOD-SUB', found); out.append(found)
        json.dump(out, open(os.path.join(HERE, 'unimod_sub_p%d.json' % p), 'w'), indent=1)
log('DONE extra3', p)
