"""(a) polished re-verification of the 5 forms with unresolved denominators; (c) best_forms.json direct check +
Hermitian-matrix reconstruction; (d) theta-series checks; (e) exact square test of zeta8_k4 idx 10."""
import numpy as np, json, sys, os, math
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vlib
HERE = os.path.dirname(os.path.abspath(__file__)); R8 = os.path.dirname(HERE)
out = {}
# (e)
m = 31965476087522569838343451499938081
print('(e) zeta8_k4 idx10 det(340*A) perfect square:', math.isqrt(m) ** 2 == m)
out['zeta8_k4_idx10_det_square'] = bool(math.isqrt(m) ** 2 == m)

# (a)
print('(a) polished verification of forms with large denominators')
out['polished'] = []
for fn in ('hur_k6_descent_s2.json', 'zeta13_k2_random_s1.json'):
    r = json.load(open(os.path.join(R8, 'data', fn)))
    name, k, n, D = r['name'], r['k'], r['n'], r['D']
    els, w, n0, extra = vlib.base_group(name)
    QB = vlib.Qbasis_reynolds(els, k)
    Qf = QB.reshape(len(QB), -1)
    gs = [np.kron(np.eye(k), g.astype(float)) for g in els]
    for i, x in enumerate(r['recs']):
        A0 = np.array(x['A_final'])
        A1 = ((Qf @ A0.ravel()) @ Qf).reshape(n, n)
        U = vlib.my_lll(A1); Uf = U.astype(float)
        Ar = Uf.T @ A1 @ Uf
        Y, _ = vlib.enum((Ar + Ar.T) / 2, 1 + 2e-3)
        X = Y @ U.T; Xf = X.astype(float)
        v = np.einsum('ij,jk,ik->i', Xf, A1, Xf)
        vs = np.sort(v)
        # gap: cluster near 1
        near = np.abs(v - 1) < 5e-5
        gap_next = vs[near.sum()] if near.sum() < len(vs) else float('nan')
        S = X[near]; Sf = S.astype(float)
        Phi = np.einsum('mi,dij,mj->md', Sf, QB, Sf, optimize=True)
        sv = np.linalg.svd(Phi, compute_uv=False)
        rank = int((sv > 1e-9 * sv[0]).sum())
        c2, *_ = np.linalg.lstsq(Phi, np.ones(len(S)), rcond=None)
        A2 = (c2 @ Qf).reshape(n, n)
        resid = np.abs(Phi @ c2 - 1).max()
        inv = max(np.abs(g.T @ A2 @ g - A2).max() for g in gs)
        U = vlib.my_lll(A2); Uf = U.astype(float)
        Ar = Uf.T @ A2 @ Uf
        Y, _ = vlib.enum((Ar + Ar.T) / 2, 1 + 1e-9)
        X2 = Y @ U.T; X2f = X2.astype(float)
        v2 = np.einsum('ij,jk,ik->i', X2f, A2, X2f)
        V2 = vlib.kappa(n) * math.exp(-0.5 * np.linalg.slogdet(A2)[1])
        d = dict(file=fn, idx=i, cluster_size=int(near.sum()), cluster_spread=float(np.abs(v[near] - 1).max()), next_value=float(gap_next),
                 rank=rank, D=D, lstsq_resid=float(resid), invariance_err=float(inv), min_polished=float(v2.min()),
                 below_1m1e9=int((v2 < 1 - 1e-9).sum()), pairs_polished=int(len(v2)), pairs_stored=x['pairs_ext'],
                 V_polished=V2, V_stored=x['V_ext'], V_relerr=abs(V2 / x['V_ext'] - 1), max_abs_diff_to_stored_A=float(np.abs(A2 - A0).max()),
                 l1_max_minvec=int(np.abs(S).sum(1).max()), ratio_best=V2 / vlib.Vbest(n))
        out['polished'].append(d)
        print('  ', json.dumps(d))

# (c) best_forms.json
print('(c) best_forms.json')
bf = json.load(open(os.path.join(R8, 'best_forms.json')))
out['best_forms'] = []
def qm(p, r):
    a, b, c, d = p; e, f, g, h = r
    return np.array([a*e - b*f - c*g - d*h, a*f + b*e + c*h - d*g, a*g - b*h + c*e + d*f, a*h + b*g - c*f + d*e])
for key, b in bf.items():
    name, k = key.split('_k'); k = int(k)
    if b['qA'] is None:
        print('  ', key, 'no integer form stored (q unresolved)'); continue
    Z = np.array(b['qA'], dtype=np.int64); q = b['q']; n = len(Z)
    detZ = vlib.bareiss_det(Z)
    V = vlib.kappa(n) * math.sqrt(q ** n / detZ)
    U = vlib.my_lll(Z); Zr = U.T @ Z @ U
    Y, _ = vlib.enum(Zr.astype(float), q + 0.5)
    X = Y @ U.T
    nrm = np.einsum('ij,jk,ik->i', X, Z, X)
    els, w, n0, extra = vlib.base_group(name)
    inv = all(np.array_equal(np.kron(np.eye(k, dtype=np.int64), g).T @ Z @ np.kron(np.eye(k, dtype=np.int64), g), Z) for g in els)
    d = dict(key=key, n=n, q=q, detZ=str(detZ), det_stated=b['det_qA'], V=V, V_stated=b['V'], min=int(nrm.min()), kissing=int(2 * (nrm == q).sum()),
             kissing_stated=b['kissing'], invariant=bool(inv), ratio_best=V / vlib.Vbest(n))
    ok = (abs(V / b['V'] - 1) < 1e-9 and d['min'] == q and d['kissing'] == b['kissing'] and inv and abs(detZ - b['det_qA']) < 1e-6 * max(1, abs(detZ)) and V <= vlib.Vbest(n) * (1 + 1e-9))
    # Hermitian reconstruction
    if 'H_times_q' in b:
        H = np.array(b['H_times_q'])
        if name == 'hur':
            beta = [np.array(v, float) for v in ([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [.5, .5, .5, .5])]
            cj = lambda p: np.array([p[0], -p[1], -p[2], -p[3]]); mul = qm
        else:
            beta = [np.array([1.0, 0.0]), np.array([-0.5, math.sqrt(3) / 2]) if name == 'eis' else np.array([0.0, 1.0])]
            cj = lambda p: np.array([p[0], -p[1]])
            mul = lambda p, r: np.array([p[0]*r[0] - p[1]*r[1], p[0]*r[1] + p[1]*r[0]])
        dd = len(beta)
        A2 = np.zeros((n, n))
        for i in range(k):
            for j in range(k):
                for a in range(dd):
                    for bb in range(dd):
                        A2[i*dd + a, j*dd + bb] = mul(mul(cj(beta[a]), beta[bb]), H[i, j])[0]
        d['H_reconstruct_err'] = float(np.abs(A2 - Z).max())
        d['H_hermitian'] = bool(all(np.abs(H[i, j] - cj(H[j, i])).max() < 1e-5 for i in range(k) for j in range(k)))
        d['H_diag'] = [float(H[i, i][0]) for i in range(k)]
        d['H_offdiag_absmax'] = float(max([np.abs(H[i, j]).max() for i in range(k) for j in range(k) if i != j] or [0]))
        ok = ok and d['H_reconstruct_err'] < 1e-4 and d['H_hermitian']
    d['ok'] = bool(ok)
    out['best_forms'].append(d)
    print('  ', ('OK  ' if ok else 'BAD ') + json.dumps(d))

# (d) theta series checks
print('(d) theta-series checks (pairs with norm 4 / 6 of the integral form 4A)')
out['theta'] = []
for key, expect in (('zeta7_k2', (378, 2016)), ('eis_k6', (378, 2016)), ('gauss_k6', (378, 2016)), ('zeta12_k3', (378, 2016)),
                    ('ico_k2', (2160, 30720)), ('hur_k4', (2160, 30720)), ('zeta8_k4', (2160, 30720)), ('zeta12_k4', (2160, 30720)),
                    ('hur_k3', None), ('zeta8_k3', None), ('hur_k5', None), ('zeta5_k3', None), ('gauss_k5', None), ('eis_k5', None)):
    b = bf[key]; Z = np.array(b['qA'], dtype=np.int64); n = len(Z)
    U = vlib.my_lll(Z); Zr = U.T @ Z @ U
    Y, _ = vlib.enum(Zr.astype(float), 6.5, cap=400000)
    X = Y @ U.T
    nrm = np.einsum('ij,jk,ik->i', X, Z, X)
    cnt = {int(t): int((nrm == t).sum()) for t in np.unique(nrm)}
    d = dict(key=key, n=n, pairs_by_norm=cnt, expected=expect, ok=(expect is None or (cnt.get(4), cnt.get(6)) == expect))
    out['theta'].append(d)
    print('  ', json.dumps(d))
json.dump(out, open(os.path.join(HERE, 'verify2.json'), 'w'), indent=1)
