"""Direct certificate for a Delta-line: build M = O(Lambda) w  (O = commutant order realised via the H-fixed lattice)
and compute R_U lower bound  r^8 D^4 / det_B(M).  No maximality assumption is needed: M is a sublattice of Lambda cap Delta w."""
import numpy as np, math, sys, json
import proto as P, run37 as R
V = P.V
def unrealify(X, p):
    h = (p - 1) // 2
    qr = np.array([P.legendre(a, p) == 1 for a in range(1, h + 1)])
    nq = int(qr.sum()); n = X.shape[0]
    Xp = X[:, :nq * 2 * h].reshape(n, nq, 2 * h); Xm = X[:, nq * 2 * h:].reshape(n, h - nq, 2 * h)
    out = np.zeros((n, h, h), dtype=complex)
    out[:, qr, :] = (Xp[:, :, :h] + 1j * Xp[:, :, h:]) / math.sqrt(2)
    out[:, ~qr, :] = (Xm[:, :, :h] + 1j * Xm[:, :, h:]) / math.sqrt(2)
    return out
def orbit_multi(S, T, vs, maxn=300000):
    Tinv = np.conj(np.transpose(T, (0, 2, 1))); gens = [S, T, Tinv]
    key = lambda u: np.round(np.concatenate([u.real.ravel(), u.imag.ravel()]) * 1e5).astype(np.int64).tobytes()
    seen = {key(vs[0]): 0}; out = [vs]; frontier = [vs]
    while frontier and len(out) < maxn:
        new = []
        for st in frontier:
            for g in gens:
                w = np.einsum('aij,kaj->kai', g, st)
                k = key(w[0])
                if k not in seen:
                    seen[k] = len(out); out.append(w); new.append(w)
        frontier = new
    return np.array(out)       # (norbit, 1+8, h, h)
def certify(p, mH, wreal, label=''):
    h = (p - 1) // 2; r = (p - 1) // 4
    for kappa in (1, -1):
        S, T = P.build(p, kappa); rel = P.check_relations(S, T, p)
        if rel['S^2=-I'] < 1e-9 and rel['(ST)^3=-I'] < 1e-9: break
    v0 = P.v0_vec(p); v0c = np.broadcast_to(v0.astype(complex), (h, h)).copy()
    Bas = np.load('Bas_p%d.npy' % p); dp = Bas.shape[1] // 2
    Bm = np.round(Bas @ Bas.T); sg, ld = np.linalg.slogdet(Bm); D = math.exp(ld / (p - 1))
    BasH = np.load('BasH_p%d_m%d.npy' % (p, mH))            # basis of Lambda^H (rank 8), H = stabiliser-type subgroup
    U = unrealify(BasH, p)
    st0 = np.concatenate([v0c[None], U])
    ORB = orbit_multi(S, T, st0)
    Xp, Xm = P.realify(ORB[:, 0], p); X = np.concatenate([Xp, Xm], axis=1)      # (norb, dim)
    c, res, rk, sv = np.linalg.lstsq(X.T, wreal, rcond=None)
    print(' %s: orbit %d, residual of w = sum c_g g v0 : %.2e' % (label, len(ORB), float(np.abs(X.T @ c - wreal).max())))
    M = []
    for i in range(1, 9):
        Yp, Ym = P.realify(ORB[:, i], p); Y = np.concatenate([Yp, Ym], axis=1)
        M.append(c @ Y)
    M = np.array(M)
    # the generated vectors must lie in Lambda: integer coordinates in Bas
    C = np.linalg.solve(Bm, Bas @ M.T).T
    print('   max distance of O.w generators from the lattice: %.2e' % float(np.abs(C - np.round(C)).max()))
    G = M @ M.T
    s8, ld8 = np.linalg.slogdet(G)
    rank = np.linalg.matrix_rank(G, tol=1e-6)
    logR = 8 * math.log(r) + 4 * math.log(D) - ld8
    mp = (wreal[:dp] ** 2).sum(); mm = (wreal[dp:] ** 2).sum(); theta = (p - 1) ** 2 / (4 * p)
    print('   rank of O.w = %d ; det_B(O.w) = exp(%.6f) ; predicted ((p/4) m+ m-)^4 = exp(%.6f)' % (rank, ld8, 4 * math.log(p / 4 * mp * mm)))
    print('   Gamma(w) = %.6f, theta = %.6f ; lower bound R_U >= r^8 D^4/det(O.w) = %.6f  (=(theta/Gamma)^4 = %.6f)' % (mp * mm / D, theta, math.exp(logR), (theta / (mp * mm / D)) ** 4))
    return math.exp(logR)
if __name__ == '__main__':
    p = int(sys.argv[1]); mH = int(sys.argv[2]); msub = int(sys.argv[3]); coords = json.loads(sys.argv[4])
    BasS = np.load('BasH_p%d_m%d.npy' % (p, msub))
    w = np.array(coords, dtype=float) @ BasS
    certify(p, mH, w, label='p=%d' % p)
