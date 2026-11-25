"""Data for the t-circle figure: min over the lattice of Q_t(w) / (2 sqrt(theta D)) as a function of t."""
import numpy as np, math, json, sys
import proto as P, run37 as R, scan as SC
V = P.V
def profile(Gp, Gm, p, D, npts=61):
    theta = (p - 1) ** 2 / (4 * p); crit = 2 * math.sqrt(theta * D)
    L = math.log(SC.fund_unit(p)[0]); ts = np.linspace(-L, L, npts); out = []
    for t in ts:
        Gt = math.exp(t) * Gp + math.exp(-t) * Gm
        b = 1.02
        while True:
            Xs, vals, flag = V.short_vectors(Gt, crit * b, maxc=400000, maxnodes=400_000_000)
            if len(Xs): out.append(float(vals.min() / crit)); break
            b *= 1.12
            if b > 3: out.append(float('nan')); break
    return [float(x / L) for x in ts], out
res = {}
# p = 13, full lattice
Bas = np.load('Bas_p13.npy'); dp = Bas.shape[1] // 2
res['13'] = profile(Bas[:, :dp] @ Bas[:, :dp].T, Bas[:, dp:] @ Bas[:, dp:].T, 13, 117.0)
# p = 37, C_3-fixed sublattice
BasH = np.load('BasH_p37_m3.npy'); dp = BasH.shape[1] // 2
res['37'] = profile(BasH[:, :dp] @ BasH[:, :dp].T, BasH[:, dp:] @ BasH[:, dp:].T, 37, 2997.0)
# p = 61, C_5-fixed sublattice (rank 24)
p = 61; m = 5; rng = np.random.default_rng(3); h = (p - 1) // 2
S, T = P.build(p, -1); O = SC.random_orbit_vectors(S, T, P.v0_vec(p), 160, rng)
g0 = next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in set(R.P_factors(p - 1))))
M = R.torus_matrix(p, pow(g0, (p - 1) // m, p)); N = np.eye(h); Mk = np.eye(h)
for _ in range(m - 1):
    Mk = Mk @ M; N = N + Mk
Yp, Ym = P.realify(np.einsum('ij,naj->nai', N, O), p)
BasH, rH = P.lattice_basis_exact(np.concatenate([Yp, Ym], axis=1)); dp = Yp.shape[1]
res['61'] = profile(BasH[:, :dp] @ BasH[:, :dp].T, BasH[:, dp:] @ BasH[:, dp:].T, 61, 15.0 ** 2 * 61)
json.dump(res, open('figdata_tcircle.json', 'w'))
for k, (ts, v) in res.items(): print(k, 'min ratio %.4f at t/L=%.2f ; max %.3f' % (np.nanmin(v), ts[int(np.nanargmin(v))], np.nanmax(v)))
