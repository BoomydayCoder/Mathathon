"""Root system of the unimodular lattice L = (C_m-fixed lattice)/sqrt(m) inside the SL_2(p) Weil lattice."""
import numpy as np, math, sys, collections
import proto as P, run37 as R, scan as SC
V = P.V
def rootsys(p, m, nsamp=160, seed=1):
    rng = np.random.default_rng(seed); h = (p - 1) // 2; r = (p - 1) // 4; D = r * r * p
    S, T = P.build(p, -1); v0 = P.v0_vec(p)
    O = SC.random_orbit_vectors(S, T, v0, nsamp, rng); del S
    g0 = next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in set(R.P_factors(p - 1))))
    a = pow(g0, (p - 1) // m, p); M = R.torus_matrix(p, a)
    N = np.eye(h); Mk = np.eye(h)
    for _ in range(m - 1):
        Mk = Mk @ M; N = N + Mk
    ON = np.einsum('ij,naj->nai', N, O)
    Yp, Ym = P.realify(ON, p)
    BasH, rH = P.lattice_basis_exact(np.concatenate([Yp, Ym], axis=1))
    dp = Yp.shape[1]
    GpH = BasH[:, :dp] @ BasH[:, :dp].T; GmH = BasH[:, dp:] @ BasH[:, dp:].T
    eps = SC.fund_unit(p)[0]; t = -math.log(eps)
    Gu = np.round((math.exp(t) * GpH + math.exp(-t) * GmH) / math.sqrt(D) / m)
    Xs, vals, flag = V.short_vectors(Gu, 2.01, maxc=100000)
    Rts = np.concatenate([Xs, -Xs]); IP = np.round(Rts @ Gu @ Rts.T).astype(int)
    n = len(Rts); seen = np.zeros(n, bool); comps = []
    for i in range(n):
        if seen[i]: continue
        st = [i]; seen[i] = True; c = [i]
        while st:
            x = st.pop()
            for y in np.where((IP[x] != 0) & (~seen))[0]:
                seen[y] = True; st.append(y); c.append(y)
        comps.append(c)
    desc = sorted((len(c), int(np.linalg.matrix_rank(IP[np.ix_(c, c)].astype(float)))) for c in comps)
    print('p=%d, C_%d-fixed lattice / sqrt(%d): rank %d, det %.3f, #roots %d, components (roots, rank): %s' % (p, m, m, rH, np.linalg.det(Gu), n, desc), flush=True)
if __name__ == '__main__':
    for a in sys.argv[1:]:
        p, m = a.split(':'); rootsys(int(p), int(m))
