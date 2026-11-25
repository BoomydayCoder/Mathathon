"""independent check of Q-extremality (= Q-perfect + strictly Q-eutactic) of the stored final forms"""
import numpy as np, json, sys, os, glob, time
from scipy.optimize import linprog
sys.path.insert(0, '.')
import vlib
t00 = time.time()
tot = 0; agree = 0; eut = 0; res = []; lpv = []
for fn in sorted(glob.glob('../data/*.json')):
    r = json.load(open(fn)); name, k, n, D = r['name'], r['k'], r['n'], r['D']
    if time.time() - t00 > 150: print('time cap, skipped', fn); continue
    els, w, n0, extra = vlib.base_group(name)
    QB = vlib.Qbasis_reynolds(els, k); Qf = QB.reshape(len(QB), -1)
    for i, x in enumerate(r['recs']):
        if 'lp_val' not in x: continue
        lpv.append(x['lp_val'])
        A = np.array(x['A_final']); q = vlib.find_q(A)
        if not q: continue
        Z = np.rint(q * A).astype(np.int64)
        U = vlib.my_lll(Z); Zr = U.T @ Z @ U
        Y, _ = vlib.enum(Zr.astype(float), q + 0.5); X = Y @ U.T
        nrm = np.einsum('ij,jk,ik->i', X, Z, X); S = X[nrm == q].astype(float)
        Phi = np.einsum('mi,dij,mj->md', S, QB, S, optimize=True)
        _, ix = np.unique(np.round(Phi, 6), axis=0, return_index=True); P = Phi[np.sort(ix)]
        g = Qf @ np.linalg.inv(Z / q).ravel()
        m = len(P)
        c = np.zeros(m + 1); c[-1] = -1
        Aeq = np.hstack([P.T, np.zeros((D, 1))])
        Aub = np.hstack([-np.eye(m), np.ones((m, 1))])
        rr = linprog(c, A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=g, bounds=[(None, None)] * m + [(None, 1.0)], method='highs')
        delta = rr.x[-1] if rr.status == 0 else float('nan')
        mine = bool(rr.status == 0 and delta > 1e-9)
        theirs = bool(x['lp_val'] > -1e-9)
        tot += 1; agree += (mine == theirs); eut += mine
        if mine != theirs: print('DISAGREE', os.path.basename(fn), i, 'delta', delta, 'status', rr.status, 'lp_val', x['lp_val'], 'orbrows', m)
        res.append(dict(file=os.path.basename(fn), idx=i, delta=float(delta), lp_val=x['lp_val'], eutactic_mine=mine))
lpv = np.array(lpv)
print('forms checked', tot, 'strictly Q-eutactic (mine)', eut, 'agree with stored LP verdict', agree, 'time %.0fs' % (time.time() - t00))
print('stored lp_val: count', len(lpv), 'negative (< -1e-9):', int((lpv < -1e-9).sum()), 'borderline |val|<1e-7:', int((np.abs(lpv) < 1e-7).sum()), 'min positive %.3e' % lpv[lpv > 1e-7].min())
json.dump(res, open('eutaxy.json', 'w'))
