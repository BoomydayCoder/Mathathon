import numpy as np, sys, math, time, json
import proto as P
V = P.V
EPS = {5: (1 + 5 ** .5) / 2, 13: (3 + 13 ** .5) / 2, 29: (5 + 29 ** .5) / 2, 37: 6 + 37 ** .5, 61: (39 + 5 * 61 ** .5) / 2}

def torus_matrix(p, a):
    """(d(a) f)(x) = (a/p) f(a x) on odd functions, basis e_x = delta_x - delta_{-x}; returns h x h real matrix."""
    h = (p - 1) // 2; M = np.zeros((h, h))
    s = P.legendre(a, p)
    # d(a) e_x = (a/p) e_{x/a}   (since (d f)(y) = f(a y), delta_x(a y)=1 iff y = x/a)
    ainv = pow(a, -1, p)
    for x in range(1, h + 1):
        y = (x * ainv) % p
        if y <= h: M[y - 1, x - 1] += s
        else: M[p - y - 1, x - 1] -= s
    return M

def analyze(p, sub_orders=(), full_lll=True, tstep=0.1, out=None):
    t0 = time.time()
    h = (p - 1) // 2; n = 2 * (p - 1)
    for kappa in (1, -1):
        S, T = P.build(p, kappa)
        rel = P.check_relations(S, T, p)
        if rel['S^2=-I'] < 1e-9 and rel['(ST)^3=-I'] < 1e-9: break
    v = P.v0_vec(p)
    O = P.orbit(S, T, v)
    print('p', p, 'kappa', kappa, 'orbit', len(O), 'time %.1f' % (time.time() - t0), flush=True)
    Xp, Xm = P.realify(O, p)
    # structured initial basis: Z[T]-span of two orbit vectors (well conditioned)
    Td = np.array([np.diag(T[a]) for a in range(h)])            # (h emb, h coords)
    def zt_span(u):
        return np.array([u * Td ** j for j in range(p - 1)])
    best = None
    rng = np.random.default_rng(1)
    cands = [O[0]] + [O[i] for i in rng.choice(len(O), size=min(12, len(O) - 1), replace=False)]
    for i1 in range(len(cands)):
        for i2 in range(i1 + 1, len(cands)):
            W = np.concatenate([zt_span(cands[i1]), zt_span(cands[i2])])
            Wp, Wm = P.realify(W, p); Wr = np.concatenate([Wp, Wm], axis=1)
            G = Wr @ Wr.T
            sg, ld = np.linalg.slogdet(G)
            if sg > 0 and np.linalg.cond(G) < 1e10 and (best is None or ld < best[0]): best = (ld, Wr)
    print(' structured initial basis: logdet %.4f' % best[0], flush=True)
    Bas, rank = P.lattice_basis(Xp, Xm, verbose=True, init=best[1])
    dp = Xp.shape[1]
    Gp = Bas[:, :dp] @ Bas[:, :dp].T; Gm = Bas[:, dp:] @ Bas[:, dp:].T
    B = Gp + Gm; E = math.sqrt(p) * (Gp - Gm)
    print(' rank', rank, '(expected %d)' % n, 'integrality', float(np.abs(B - np.round(B)).max()), float(np.abs(E - np.round(E)).max()), flush=True)
    Bi = np.round(B); Ei = np.round(E)
    print(' E B^-1 E - pB residual', float(np.abs(Ei @ np.linalg.solve(Bi, Ei) - p * Bi).max()))
    sign, logdet = np.linalg.slogdet(Bi)
    D = math.exp(logdet / (p - 1)); theta = (p - 1) ** 2 / (4 * p)
    print(' logdet B = %.6f, D = %.8f, theta = %.6f, critical Q-bound 2 sqrt(theta D) = %.6f' % (logdet, D, theta, 2 * math.sqrt(theta * D)), flush=True)
    np.save('B_p%d.npy' % p, Bi); np.save('E_p%d.npy' % p, Ei); np.save('Bas_p%d.npy' % p, Bas)
    L = math.log(EPS[p]); ts = np.arange(-L, L + 1e-9, tstep)
    bound = 2 * math.sqrt(theta * D) * math.cosh(tstep / 2) * (1 + 1e-9)
    res = {'p': p, 'rank': rank, 'logdet': logdet, 'D': D, 'theta': theta, 'sub': {}}
    # --- fixed sublattices of odd-order split torus subgroups
    X = np.concatenate([Xp, Xm], axis=1)
    g0 = None
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in set(P_factors(p - 1))): g0 = g; break
    for m in sub_orders:
        a = pow(g0, (p - 1) // m, p)
        M = torus_matrix(p, a)
        # norm map on the complex orbit vectors (same matrix in every embedding)
        N = np.eye(h); Mk = np.eye(h)
        for _ in range(m - 1):
            Mk = Mk @ M; N = N + Mk
        ON = np.einsum('ij,naj->nai', N, O)
        Yp, Ym = P.realify(ON, p)
        keep = np.where((np.abs(Yp).sum(1) + np.abs(Ym).sum(1)) > 1e-9)[0]
        Yp, Ym = Yp[keep], Ym[keep]
        try:
            BasH, rH = P.lattice_basis(Yp, Ym, verbose=False)
        except Exception as e:
            print(' sub order', m, 'lattice_basis failed', e); continue
        GpH = BasH[:, :dp] @ BasH[:, :dp].T; GmH = BasH[:, dp:] @ BasH[:, dp:].T
        s2, ld = np.linalg.slogdet(np.round(GpH + GmH))
        r = (p - 1) // 4
        # R_U = (n/rH)^{rH} det(Lambda)^{rH/n} / det(Lambda_U)   (norm-map image; a sublattice of the fixed lattice)
        logR = rH * math.log(n / rH) + logdet * rH / n - ld
        print(' subgroup C_%d: rank of N_H(Lambda) = %d, logdet = %.4f, log R_U(image of norm map) = %.6f' % (m, rH, ld, logR), flush=True)
        best = (1e99, None); tot = 0; t1 = time.time()
        for t in ts:
            Gt = math.exp(t) * GpH + math.exp(-t) * GmH
            try:
                Xs, vals, flag = V.short_vectors(Gt, bound, maxc=400000, maxnodes=400_000_000)
            except RuntimeError as e:
                print('   t=%.2f node limit' % t); continue
            if len(Xs):
                mp = np.einsum('ij,jk,ik->i', Xs, GpH, Xs); mm = np.einsum('ij,jk,ik->i', Xs, GmH, Xs)
                gam = mp * mm / D; tot += len(Xs)
                if gam.min() < best[0]: best = (float(gam.min()), float(t))
        print('   exhaustive over t-grid: #vectors below bound (sum over t) = %d, min Gamma = %s at t=%s  (theta = %.6f)  [%.1fs]' % (tot, best[0], best[1], theta, time.time() - t1), flush=True)
        res['sub'][m] = {'rank': rH, 'logR': logR, 'gamma_min': best[0], 'count': tot}
    # --- heuristic: LLL on the full lattice over the t-grid
    if full_lll:
        best = (1e99, None, None); bestQ = (1e99, None)
        for t in ts:
            Gt = math.exp(t) * Gp + math.exp(-t) * Gm
            G2, U = V.lll_gram(Gt)
            Xs = U.T
            mp = np.einsum('ij,jk,ik->i', Xs, Gp, Xs); mm = np.einsum('ij,jk,ik->i', Xs, Gm, Xs)
            gam = mp * mm / D
            i = int(np.argmin(gam))
            if gam[i] < best[0]: best = (float(gam[i]), float(t), Xs[i].copy())
            qn = np.diag(G2).min() / (2 * math.sqrt(theta * D))
            if qn < bestQ[0]: bestQ = (float(qn), float(t))
        print(' full lattice, LLL bases over t-grid: min Gamma among basis vectors = %.6f (theta %.6f, ratio %.4f) at t=%.2f ; min Q_t/critical = %.4f' % (best[0], theta, best[0] / theta, best[1], bestQ[0]), flush=True)
        res['lll_gamma_min'] = best[0]; res['lll_ratio'] = best[0] / theta
        if best[0] < theta * (1 - 1e-9):
            json.dump({'p': p, 'w': [int(x) for x in best[2]], 'gamma': best[0], 'theta': theta, 'D': D, 't': best[1]}, open('violator_p%d.json' % p, 'w'))
            print(' *** VIOLATOR FOUND (certificate saved) ***')
    print(' total time %.1f' % (time.time() - t0), flush=True)
    if out: json.dump(res, open(out, 'w'))
    return res

def P_factors(n):
    f = []; d = 2
    while d * d <= n:
        while n % d == 0: f.append(d); n //= d
        d += 1
    if n > 1: f.append(n)
    return f

if __name__ == '__main__':
    p = int(sys.argv[1]); subs = tuple(int(x) for x in sys.argv[2].split(',')) if len(sys.argv) > 2 and sys.argv[2] else ()
    analyze(p, subs, full_lll=(len(sys.argv) <= 3 or sys.argv[3] != 'nolll'), out='res_p%d.json' % p)
