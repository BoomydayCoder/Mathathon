#!/usr/bin/env python3
# Independent end-to-end verification of saved lattice data and violator certificates (Builder B).
import sys, json, math, os
import numpy as np
from weil import *
p = int(sys.argv[1]); certs = sys.argv[2:]
Wl = Weil(p); n, m, h, Dn = Wl.n, Wl.m, Wl.h, Wl.Dn
B = np.load(os.path.join(HERE, 'B_p%d.npy' % p)); E = np.load(os.path.join(HERE, 'E_p%d.npy' % p)); Ar = np.load(os.path.join(HERE, 'Ared_p%d.npy' % p))
out = dict(p=p)
# (1) Gram matrices recomputed numerically from the complex embeddings of the saved basis
V = (Ar.astype(float) / Dn) @ Wl.Bas.reshape(n, -1)            # n x (h*h) complex ; rows = sigma_a(u)_x
V3 = V.reshape(n, h, h)
plus = Wl.chi == 1
Gp = 2 * np.real(V3[:, plus, :].reshape(n, -1).conj() @ V3[:, plus, :].reshape(n, -1).T)
Gm = 2 * np.real(V3[:, ~plus, :].reshape(n, -1).conj() @ V3[:, ~plus, :].reshape(n, -1).T)
out['max|Gp+Gm-B|'] = float(np.abs(Gp + Gm - B).max()); out['max|sqrt(p)(Gp-Gm)-E|'] = float(np.abs(math.sqrt(p) * (Gp - Gm) - E).max())
log('numeric Gram check: |G+ + G- - B| = %.2e, |sqrt p (G+ - G-) - E| = %.2e' % (out['max|Gp+Gm-B|'], out['max|sqrt(p)(Gp-Gm)-E|']))
# (2) exact determinant of saved B (Bareiss)
dB = det_exact(B); out['detB_equals_(m^2 p)^(p-1)'] = (dB == (m * m * p) ** (p - 1))
log('exact det(B_saved) == (m^2 p)^(p-1): %s' % out['detB_equals_(m^2 p)^(p-1)'])
# (3) orbit of v0 by BFS in the complex model; coordinates in the saved basis; Z-span = whole lattice ?
v0 = np.tile(Wl.Wc[0][None, :], (h, 1)).astype(complex)       # sigma_a(v0) = v0 (rational)
def key(v): return tuple(np.round(v[0].real * 1e5).astype(np.int64)) + tuple(np.round(v[0].imag * 1e5).astype(np.int64))
seen = {key(v0): 0}; orbit = [v0]; front = [v0]
while front:
    nf = []
    for v in front:
        for w in (np.einsum('axy,ay->ax', Wl.S_all, v), Wl.Td * v):
            kk = key(w)
            if kk not in seen: seen[kk] = len(orbit); orbit.append(w); nf.append(w)
    front = nf
    if len(orbit) > p * (p * p - 1): break
out['orbit_size'] = len(orbit); out['expected_orbit_size_|G|/m'] = p * (p * p - 1) // m
Ov = np.array(orbit).reshape(len(orbit), -1)
Y = 2 * np.real(Ov.conj() @ V.T)                               # B(orbit vec, basis)
Cf = np.linalg.solve(B.astype(float), Y.T).T; Ci = np.rint(Cf); res = float(np.abs(Cf - Ci).max())
Hn = hnf_mod(Ci.astype(np.int64), Dn)
out['orbit_coords_integrality_residual'] = res; out['orbit_spans_saved_lattice'] = bool(np.all(np.diag(Hn) == 1))
log('orbit of v0: size %d (expected |G|/|H| = %d); coordinates in saved basis integral (residual %.2e); Z-span of orbit == saved lattice: %s' %
    (len(orbit), out['expected_orbit_size_|G|/m'], res, out['orbit_spans_saved_lattice']))
# (4) certificates
for cf in certs:
    c = json.load(open(os.path.join(HERE, cf))); w = np.array(c['w_coords_in_saved_basis'], dtype=object)
    b = int(w @ B.astype(object) @ w); e = int(w @ E.astype(object) @ w); N = p * b * b - e * e; thr = m * m * p * (p - 1) ** 2
    wv = (np.array(c['w_coords_in_saved_basis'], dtype=float) @ V).reshape(h, h)
    mp_ = float(2 * (np.abs(wv[plus]) ** 2).sum()); mm_ = float(2 * (np.abs(wv[~plus]) ** 2).sum())
    D = m * m * p; gam = mp_ * mm_ / D
    r = dict(file=cf, b=b, e=e, N=N, threshold=thr, exact_violator=(N < thr), m_plus_direct=mp_, m_minus_direct=mm_, Gamma_direct=gam, Gamma_exact='%d/%d' % (N, 4 * D * p),
             theta=Wl.theta, ratio=gam / Wl.theta, R_U_lower_bound=(Wl.theta / gam) ** 4)
    log('CERT CHECK', r); out[cf] = r
json.dump(out, open(os.path.join(HERE, 'verify_p%d.json' % p), 'w'), indent=1)
