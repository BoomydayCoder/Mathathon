#!/usr/bin/env python3
"""Independent re-verification of every saved form in the R9 data (own enumeration code, own LLL).
Usage: python3 verify_all.py out.json file1.json file2.json ..."""
import sys, os, json, math, time
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vlib import *

def eutaxy_margin(A, S):
    n = A.shape[0]; iu = np.triu_indices(n)
    Phi = (S[:, iu[0]] * S[:, iu[1]]).astype(float); m = len(S)
    Aeq = np.hstack([Phi.T, np.zeros((len(iu[0]), 1))]); beq = np.linalg.inv(A)[iu]
    Aub = np.hstack([-np.eye(m), np.ones((m, 1))]); c = np.zeros(m + 1); c[-1] = -1
    res = linprog(c, A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=beq, bounds=[(None, None)] * (m + 1), method='highs')
    return float(res.x[-1]) if res.status == 0 else float('nan')

out = sys.argv[1]; files = sys.argv[2:]
res = []; t00 = time.time()
for f in files:
    d = json.load(open(f)); n = d['n']; N = n * (n + 1) // 2; kind = d['kind']
    Vb = Vbest(n)
    t0 = time.time(); nforms = 0
    for r in d['recs']:
        o = dict(file=os.path.basename(f), n=n, kind=kind, idx=r['idx'], has_lp='A_LP' in r, errors=[k for k in r if 'error' in k])
        Q = np.array(r['Q'])
        o['Q_sym'] = float(np.abs(Q - Q.T).max()); o['Q_det'] = float(np.linalg.det(Q)); o['Q_mineig'] = float(np.linalg.eigvalsh(Q)[0])
        X, v, _ = short(Q, Q.diagonal().min() * (1 + 1e-12))
        lam2 = float(v.min()); o['lam1sq'] = lam2; o['V0'] = kap(n) * lam2 ** (n / 2)
        o['lam1sq_stored'] = r.get('lam1sq'); o['V0_stored'] = r.get('V0')
        Qi = np.linalg.inv(Q)
        for tag in ('LP', 'K', 'Kd', 'extLP', 'extK'):
            if 'A_' + tag not in r: continue
            A = np.array(r['A_' + tag]); nforms += 1
            mn, npairs, S = min_and_pairs(A)
            V = kap(n) / math.sqrt(np.linalg.det(A))
            e = dict(min=mn, pairs=npairs, V=V, V_stored=r['V_' + tag], pairs_stored=r['pairs_' + tag], rank=rank_perfection(S),
                     V_over_Vbest=V / Vb, sym=float(np.abs(A - A.T).max()))
            if tag in ('LP', 'K', 'Kd'):
                w = np.linalg.eigvals(Qi @ A).real; m = np.trace(Qi @ A) / n
                e['v'] = kap(n) * m ** (-n / 2); e['v_stored'] = r['v_' + tag]
                e['B0sq'] = float(((np.sort(w) / m - 1) ** 2).sum()); e['B0sq_stored'] = r['B0sq' if tag == 'LP' else 'B0sq_' + tag]
                e['V_viaQ'] = kap(n) * math.sqrt(np.linalg.det(Q) / np.linalg.det(A))
            if tag == 'LP':
                e['eig_ratio'] = float(w.max() / w.min()); e['eig_ratio_stored'] = r['eig_ratio']
                rr, sl, npos, lam = dual_certificate(Q, S)
                e['dual_resid'] = float(rr); e['dual_sum'] = sl; e['primal_obj'] = float(np.trace(Qi @ A)); e['lp_obj_stored'] = r['lp_info']['lp_obj']
                e['dual_npos'] = npos; e['dual_minpos'] = float(lam[lam > 1e-9].min())
                e['rounds'] = r['lp_info']['rounds']; e['fallback'] = r['lp_info']['fallback']; e['enlarge'] = r['lp_info']['enlarge']; e['ncons'] = r['lp_info']['ncons']
                e['rank_stored'] = r['rank_LP']
            if tag in ('extLP', 'extK'):
                e['eutaxy'] = eutaxy_margin(A, S); e['steps'] = r['steps_' + tag]
            if tag in ('K', 'Kd'):
                e['steps'] = r['steps_' + tag]
            o[tag] = e
        res.append(o)
    print(f'{os.path.basename(f)}: {len(d["recs"])} recs, {nforms} forms, {time.time()-t0:.1f}s', flush=True)
    json.dump(res, open(out, 'w'))
print('total', time.time() - t00)
