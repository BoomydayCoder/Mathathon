#!/usr/bin/env python3
# is the member Q_t/sqrt(D) at t = +-log(eps) integral / even / unimodular ?  (uses saved exact B, E)
import sys, json, math, os
import numpy as np
from weil import *
p = int(sys.argv[1]); bound = float(sys.argv[2]); maxnodes = int(float(sys.argv[3]))
B = np.load(os.path.join(HERE, 'B_p%d.npy' % p)); E = np.load(os.path.join(HERE, 'E_p%d.npy' % p)); m = (p - 1) // 4; D = m * m * p; n = len(B)
eps, (x, y, s) = fund_unit(p); sp = math.sqrt(p); res = dict(p=p, eps_xy_sign=(x, y, s))
# exact: Q = (eps G+ + eps^-1 G-)/(m sqrt p), G+- = (B +- E/sqrt p)/2, eps = (x + y sqrt p)/2, eps^-1 = -s/4*(x - y sqrt p)/2 *... use eps*eps' = s/4*4 -> eps' = (x - y sqrt p)/2, eps*eps' = (x^2 - p y^2)/4 = s/4
# for s = -4: eps^-1 = -eps'. Then eps G+ + eps^-1 G- = eps G+ - eps' G- = [ (x+y sp)(B + E/sp) - (x - y sp)(B - E/sp) ]/4 = [ 2 y sp B + 2 x E/sp ]/4 ; divide by m sp: (y p B + x E)/(2 m p)
if s == -4:
    num = y * p * B.astype(object) + x * E.astype(object); den = 2 * m * p
    integral = not np.any(num % den); res['Q_exact_integral'] = bool(integral)
    if integral:
        Q = (num // den); res['even'] = bool(all(int(Q[i, i]) % 2 == 0 for i in range(n))); res['diag_min'] = int(min(Q[i, i] for i in range(n)))
        Qi = np.array(Q, dtype=np.int64); res['symmetric'] = bool(np.array_equal(Qi, Qi.T))
        sgn, ld = np.linalg.slogdet(Qi.astype(float)); res['logdet_float'] = float(ld)
        if n <= 80: res['det_exact'] = str(det_exact(Qi))
        U, Qr = lll_int(Qi, 0.999); res['lll_diag_min'] = int(Qr.diagonal().min())
        if bound > 0:
            X, vals, flag, nodes = fp(Qr.astype(float), bound, maxc=1000000, maxnodes=maxnodes)
            res['enum'] = dict(bound=bound, pairs_found=len(X), flag=flag, nodes=nodes, minval=float(vals.min()) if len(vals) else None)
log('UNIMOD-EXACT', res)
json.dump(res, open(os.path.join(HERE, 'unimod_p%d.json' % p), 'w'), indent=1)
