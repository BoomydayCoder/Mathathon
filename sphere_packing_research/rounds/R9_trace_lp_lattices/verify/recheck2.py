"""Second pass: (1) explicit dual-certificate residuals for all LP vertices; (2) all walk/greedy forms re-examined in
LLL-reduced coordinates of A itself after re-polishing (to separate storage rounding from genuine failures)."""
import json, numpy as np, math, sys, glob, os
sys.path.insert(0, '.')
from vlib import *
from scipy.optimize import linprog, nnls
def eutaxy_margin(A, S):
    n = A.shape[0]; iu = np.triu_indices(n)
    Phi = (S[:, iu[0]] * S[:, iu[1]]).astype(float); m = len(S)
    Aeq = np.hstack([Phi.T, np.zeros((len(iu[0]), 1))]); beq = np.linalg.inv(A)[iu]
    Aub = np.hstack([-np.eye(m), np.ones((m, 1))]); c = np.zeros(m + 1); c[-1] = -1
    res = linprog(c, A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=beq, bounds=[(None, None)] * (m + 1), method='highs')
    return float(res.x[-1]) if res.status == 0 else float('nan')
files = sorted(glob.glob('../data/n*_*.json')) + sorted(glob.glob('../data/lpscan_*.json'))
worst = dict(resid=0, gap=0, neg=0); nlp = 0; minlam = 1e9; ndeg = 0
out = []
for f in files:
    d = json.load(open(f)); n = d['n']; N = n * (n + 1) // 2; iu = np.triu_indices(n)
    for r in d['recs']:
        Q = np.array(r['Q']); Qi = np.linalg.inv(Q)
        A = np.array(r['A_LP']); mn, npairs, S = min_and_pairs(A)
        Phi = (S[:, iu[0]] * S[:, iu[1]]).astype(float); c = Qi[iu]
        lam, _ = nnls(Phi.T, c, maxiter=100000)
        resid = np.linalg.norm(Phi.T @ lam - c) / np.linalg.norm(c)
        gap = abs(lam.sum() - np.trace(Qi @ A)) / np.trace(Qi @ A)
        nlp += 1; worst['resid'] = max(worst['resid'], resid); worst['gap'] = max(worst['gap'], gap)
        pos = lam[lam > 1e-12]; minlam = min(minlam, pos.min()); ndeg += (len(pos) < N)
        # exact-arithmetic style cross-check for simple vertices: unique multipliers by linear solve
        if npairs == N:
            lam2 = np.linalg.solve(Phi.T, c); worst['neg'] = min(worst['neg'], lam2.min()) if lam2.min() < 0 else worst['neg']
        for tag in ('K', 'Kd', 'extLP', 'extK'):
            if 'A_' + tag not in r: continue
            A = np.array(r['A_' + tag])
            G, U = my_lll(A)                     # reduced coordinates of the form itself
            X, v, _ = short(G, 1 + 1e-6)
            M = np.einsum('ki,kj->kij', X, X).reshape(len(X), n * n).astype(float)
            dlt = np.linalg.lstsq(M, 1.0 - M @ G.ravel(), rcond=None)[0].reshape(n, n); dlt = (dlt + dlt.T) / 2
            G2 = G + dlt
            X2, v2, _ = short(G2, 1 + 1e-6)
            S2 = X2[np.abs(v2 - 1) < 1e-9]
            Ui = np.linalg.inv(U.astype(float))
            dA = Ui.T @ dlt @ Ui                 # correction expressed in the stored coordinates
            V2 = kap(n) / math.sqrt(np.linalg.det(G2))
            o = dict(file=os.path.basename(f), idx=r['idx'], tag=tag, n=n, min=float(v2.min()), pairs=len(S2), pairs_stored=r['pairs_' + tag], rank=rank_perfection(S2),
                     corr_reduced=float(np.abs(dlt).max()), corr_stored_coords=float(np.abs(dA).max()), V=V2, V_stored=r['V_' + tag], VoverVbest=V2 / Vbest(n),
                     maxU=int(np.abs(U).max()), steps=r.get('steps_' + tag))
            if tag.startswith('ext'): o['eutaxy'] = eutaxy_margin(G2, S2)
            out.append(o)
print('LP vertices: %d ; max explicit dual residual ||sum lam xx^T - Q^-1||/||Q^-1|| = %.2e ; max duality gap = %.2e ; min positive multiplier = %.2e ; #with fewer than N positive multipliers = %d ; most negative multiplier among simple vertices (linear solve) = %.2e' % (nlp, worst['resid'], worst['gap'], minlam, ndeg, worst['neg']))
json.dump(out, open('recheck2_forms.json', 'w'))
for tag in ('K', 'Kd', 'extLP', 'extK'):
    E = [o for o in out if o['tag'] == tag]
    print(f'--- {tag} ({len(E)} forms, re-polished in own reduced coordinates): max |min-1| = %.2e ; pairs mismatches %d ; rank<N: %d ; max correction (stored coords) = %.2e ; (reduced coords) = %.2e ; max |V/V_stored-1| = %.2e ; max V/Vbest = %.10f ; max |U| = %d' % (
        max(abs(o['min'] - 1) for o in E), sum(o['pairs'] != o['pairs_stored'] for o in E), sum(o['rank'] != o['n'] * (o['n'] + 1) // 2 for o in E),
        max(o['corr_stored_coords'] for o in E), max(o['corr_reduced'] for o in E), max(abs(o['V'] / o['V_stored'] - 1) for o in E), max(o['VoverVbest'] for o in E), max(o['maxU'] for o in E)))
    if tag.startswith('ext'):
        bad = [(o['file'], o['idx'], round(o['eutaxy'], 6), o['steps'], o['pairs'], round(o['V'] / (o['n'] * (o['n'] + 1) / 2), 4)) for o in E if not o['eutaxy'] > 1e-9]
        print('    non-eutactic (NOT extreme):', bad, ' min positive eutaxy margin %.3e' % min(o['eutaxy'] for o in E if o['eutaxy'] > 1e-9))
