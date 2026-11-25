#!/usr/bin/env python3
import sys, json, math, os, time
import numpy as np
from weil import *
p = int(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else 'all'
Wl, out = run(p, [])
n, m = Wl.n, Wl.m
res = dict(basic=out)
def save():
    json.dump(res, open(os.path.join(HERE, 'results_p%d_%s.json' % (p, mode + ('_' + sys.argv[5] if len(sys.argv) > 5 else ''))), 'w'), indent=1, default=str)
def calib(C, label):
    Bs, Es = Wl.grams(C); k = len(C); d = det_exact(Bs)
    R = math.exp(k * math.log(n / k) + (k / 2) * math.log(Wl.D) - math.log(d))
    exact = None
    if Wl.Dexact and n % k == 0 and k % 2 == 0: exact = (d == (n // k) ** k * Wl.Dexact ** (k // 2))
    log('CALIB %s: rank %d det_B=%s R_U=%.9f exactly borderline (R_U=1): %s' % (label, k, d if d < 10**30 else '~1e%d' % (len(str(d)) - 1), R, exact))
    return dict(label=label, rank=k, det=str(d), R_U=R, exactly_borderline=exact)
def cert(r, label):
    """exact integer certificate check for the best vector of a search result r"""
    w = np.array(r['w'], dtype=np.int64)
    b = int(w @ Wl.B @ w); e = int(w @ Wl.E @ w); N = p * b * b - e * e
    sp = math.sqrt(p); mp_, mm_ = (b + e / sp) / 2, (b - e / sp) / 2
    c = dict(p=p, label=label, w_coords_in_saved_basis=w.tolist(), b=b, e=e, m_plus=mp_, m_minus=mm_, D=Wl.D, D_exact=Wl.Dexact,
             Gamma=N / (4 * Wl.D * p), theta=Wl.theta, N_exact=N, threshold_exact=(Wl.Dexact * (p - 1) ** 2 if Wl.Dexact else None),
             is_violator_exact=(N < Wl.Dexact * (p - 1) ** 2) if Wl.Dexact else None, lower_bound_R_U=(Wl.theta / (N / (4 * Wl.D * p))) ** 4)
    return c
def unimod(t, bound, label):
    sp = math.sqrt(p); Bf = Wl.B.astype(float); Ef = Wl.E.astype(float)
    Q = (math.exp(t) * (Bf + Ef / sp) / 2 + math.exp(-t) * (Bf - Ef / sp) / 2) / math.sqrt(Wl.D)
    resid = float(np.abs(Q - np.rint(Q)).max()); even = bool(np.all(np.rint(Q.diagonal()) % 2 == 0))
    U, _ = lll_int(np.rint(Q * 2.0 ** 30).astype(np.int64)); Uf = U.astype(float); Q2 = Uf @ Q @ Uf.T
    X, vals, flag, nodes = fp(Q2, bound, maxc=300000)
    r = dict(label=label, t=t, integrality_residual=resid, even_diag=even, bound=bound, pairs_found=len(X), minval=float(vals.min()) if len(vals) else None, flag=flag)
    log('UNIMOD', r); return r
le = math.log(Wl.eps)
I = np.eye(n, dtype=np.int64)
if p == 5:
    res['calib_full'] = calib(I, 'full')
    res['full'] = Wl.tsearch(I, label='p5 full exhaustive')
    res['unimod'] = [unimod(t, 2 + 1e-6, 't=%s' % s) for t, s in ((0.0, '0'), (le, '+log eps'), (-le, '-log eps'))]
if p == 13:
    res['calib_H'] = calib(Wl.CH, 'C_3 line')
    res['H'] = Wl.tsearch(Wl.CH, label='p13 C_3 line exhaustive')
    res['full'] = Wl.tsearch(I, label='p13 full rank-24 exhaustive')
    res['unimod'] = [unimod(t, 4 + 1e-6, 't=%s' % s) for t, s in ((0.0, '0'), (le, '+log eps'), (-le, '-log eps'))]
if p == 29:
    res['calib_H'] = calib(Wl.CH, 'C_7 split line'); res['H'] = Wl.tsearch(Wl.CH, label='p29 C_7 split line exhaustive')
    P5 = Wl.word('S' + 'T' * 5); Q5 = P5
    for _ in range(4): Q5 = Q5 @ P5
    assert np.array_equal(Q5, I), 'ST^5 not of order 5'
    C5 = Wl.fixed(P5); res['calib_ns5'] = calib(C5, 'non-split order 5 fixed')
    res['ns5'] = Wl.tsearch(C5, label='p29 nonsplit C_5 fixed exhaustive'); res['cert_ns5'] = cert(res['ns5'], 'nonsplit C_5 (S T^5)')
    log('CERT', {k: v for k, v in res['cert_ns5'].items() if k != 'w_coords_in_saved_basis'})
    save()
    P3 = Wl.word('STST'); assert np.array_equal(P3 @ P3 @ P3, I)
    C3 = Wl.fixed(P3); res['calib_ns3'] = calib(C3, 'non-split order 3 fixed')
    res['ns3'] = Wl.tsearch(C3, label='p29 nonsplit C_3 fixed exhaustive'); res['cert_ns3'] = cert(res['ns3'], 'nonsplit C_3 ((ST)^2)')
    log('CERT', {k: v for k, v in res['cert_ns3'].items() if k != 'w_coords_in_saved_basis'})
    if res['cert_ns3']['is_violator_exact']: json.dump(res['cert_ns3'], open(os.path.join(HERE, 'violator_p29.json'), 'w'), indent=1)
if p in (37, 61):
    qs = {37: [9, 3], 61: [15, 5, 3]}[p]
    if mode in ('all', 'small'):
        for q in qs:
            if (p, q) == (61, 3): continue
            C = Wl.CH if q == m else Wl.fixed_perm(q)
            res['calib_C%d' % q] = calib(C, 'split C_%d fixed' % q)
            res['C%d' % q] = Wl.tsearch(C, label='p%d split C_%d fixed exhaustive' % (p, q))
            res['cert_C%d' % q] = cert(res['C%d' % q], 'split C_%d' % q); save()
            if res['cert_C%d' % q]['is_violator_exact']:
                json.dump(res['cert_C%d' % q], open(os.path.join(HERE, 'violator_p%d_C%d.json' % (p, q)), 'w'), indent=1)
                log('VIOLATOR CERT', {k: v for k, v in res['cert_C%d' % q].items() if k != 'w_coords_in_saved_basis'})
    if mode == 'c3big':
        C = Wl.fixed_perm(3)
        res['calib_C3'] = calib(C, 'split C_3 fixed'); save()
        res['C3'] = Wl.tsearch(C, label='p%d split C_3 fixed (rank 40) exhaustive' % p, maxnodes=int(float(sys.argv[3])), tlim=float(sys.argv[4]), trange=(int(sys.argv[5]), int(sys.argv[6])))
        res['cert_C3'] = cert(res['C3'], 'split C_3'); save()
        log('CERT', {k: v for k, v in res['cert_C3'].items() if k != 'w_coords_in_saved_basis'})
    if mode in ('all', 'full'):
        res['full_heur'] = Wl.tsearch(I, do_enum=False, pairs=True, label='p%d FULL lattice heuristic (LLL over t-grid + pair sums)' % p)
        res['cert_full'] = cert(res['full_heur'], 'full heuristic'); save()
        log('CERT', {k: v for k, v in res['cert_full'].items() if k != 'w_coords_in_saved_basis'})
        if res['cert_full']['is_violator_exact']:
            json.dump(res['cert_full'], open(os.path.join(HERE, 'violator_p%d_full.json' % p), 'w'), indent=1)
save()
log('DONE')
