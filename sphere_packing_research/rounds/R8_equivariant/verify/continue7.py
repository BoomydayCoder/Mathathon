"""The 7 final forms whose LP value is 0 (weakly eutactic): continue the greedy descent with a strict threshold and
verify the end points with my own enumerator."""
import numpy as np, json, sys, os, math
sys.path.insert(0, '.'); sys.path.insert(0, '..')
import vlib, eq as their
def greedy2(self, c, S, rng, maxit=200):
    nst = 0; val = float('nan')
    while nst < maxit:
        val, b, r, norb = self.descent_lp(c, S)
        if val > 1e-9: break
        V0 = self.V(c)
        try: c2, S2, t = self.line_step(c, b)
        except ValueError: break
        if self.null(self.phi(S2))[0] < self.D: c2, S2, _ = self.to_perfect(c2, S2, rng, 'descent')
        if not (self.V(c2) > V0 * (1 + 1e-10)): break
        c, S = c2, S2; nst += 1
    return c, S, nst, val
out = []
for fn, i in (('gauss_k5_random_s1', 35), ('gauss_k6_random_s1', 6), ('hur_k4_random_s1', 31), ('zeta12_k3_random_s1', 15),
              ('zeta12_k3_random_s1', 17), ('zeta12_k3_random_s1', 24), ('zeta8_k3_random_s1', 24)):
    r = json.load(open(f'../data/{fn}.json')); x = r['recs'][i]; name, k, n = r['name'], r['k'], r['n']
    E = their.Eq(name, k)
    A = np.array(x['A_final']); q = vlib.find_q(A); A = np.rint(q * A) / q
    c = E.QBf @ A.ravel()
    S, v, fl = E.short(E.mat(c), 1 + 1e-8, maxc=120000)
    c = E.polish(c, E.phi(S))
    rng = np.random.default_rng(12345)
    c2, S2, nst, val = greedy2(E, c, S, rng)
    A2 = E.mat(c2)
    # my verification of the end point
    q2 = 0
    for qq in range(1, 3001):
        if np.abs(qq * A2 - np.rint(qq * A2)).max() < 1e-8 * qq: q2 = qq; break
    if q2:
        Z = np.rint(q2 * A2).astype(np.int64); U = vlib.my_lll(Z); Zr = U.T @ Z @ U
        Y, _ = vlib.enum(Zr.astype(float), q2 + 0.5); X = Y @ U.T; nrm = np.einsum('ij,jk,ik->i', X, Z, X)
        mn = int(nrm.min()) / q2; pairs = int((nrm == q2).sum()); V2 = vlib.kappa(n) * math.sqrt(q2 ** n / vlib.bareiss_det(Z))
    else:
        U = vlib.my_lll(A2); Uf = U.astype(float); Ar = Uf.T @ A2 @ Uf; Y, _ = vlib.enum((Ar + Ar.T) / 2, 1 + 1e-8)
        X = (Y @ U.T).astype(float); nrm = np.einsum('ij,jk,ik->i', X, A2, X); mn = float(nrm.min()); pairs = len(nrm)
        V2 = vlib.kappa(n) * math.exp(-.5 * np.linalg.slogdet(A2)[1])
    d = dict(file=fn, idx=i, n=n, V_stored_ext=x['V_ext'], pairs_stored=x['pairs_ext'], V_after_continuation=V2, extra_steps=nst, final_lp_val=val,
             min_check=mn, pairs=pairs, q=q2, ratio_best_before=x['V_ext'] / vlib.Vbest(n), ratio_best_after=V2 / vlib.Vbest(n))
    out.append(d); print(json.dumps(d), flush=True)
json.dump(out, open('continue7.json', 'w'), indent=1)
