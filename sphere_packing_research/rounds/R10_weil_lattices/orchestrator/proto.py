"""Orchestrator's prototype: SL_2(p) odd Weil lattice Lambda_0 = Z[G] v0 via Minkowski embedding.
Checks the specification handed to the round-10 builders (independent third implementation)."""
import numpy as np, sys, math, time
from fractions import Fraction
sys.path.insert(0, '/home/user/Mathathon/sphere_packing_research/rounds/R2_voronoi_descent')
import voronoi as V

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def build(p, kappa=1):
    h = (p - 1) // 2
    xs = np.arange(1, h + 1)
    S = np.zeros((h, h, h), dtype=complex); T = np.zeros((h, h, h), dtype=complex)
    for ia, a in enumerate(range(1, h + 1)):
        z = np.exp(2j * np.pi * a / p)
        T[ia] = np.diag(z ** ((xs * xs) % p))
        e = (2 * np.outer(xs, xs)) % p
        S[ia] = kappa * (z ** e - z ** (-e)) / (legendre(a, p) * math.sqrt(p))
    return S, T

def check_relations(S, T, p):
    h = S.shape[1]; I = np.eye(h)
    r = {}
    S2 = S @ S
    r['S^2=-I'] = float(np.abs(S2 + I).max())
    ST = S @ T
    ST3 = ST @ ST @ ST
    r['(ST)^3=-I'] = float(np.abs(ST3 + I).max()); r['(ST)^3=+I'] = float(np.abs(ST3 - I).max())
    Tp = np.linalg.matrix_power(T[0], p); r['T^p=I'] = float(np.abs(Tp - I).max())
    return r

def v0_vec(p):
    h = (p - 1) // 2
    v = np.zeros(h)
    if p == 5:
        v[0] = 1; return v
    assert p % 8 == 5
    m = (p - 1) // 4
    for x in range(1, p):
        if pow(x, m, p) == 1:
            if x <= h: v[x - 1] += 1
            else: v[p - x - 1] -= 1
    return v

def orbit(S, T, v, maxn=200000):
    Tinv = np.conj(np.transpose(T, (0, 2, 1)))
    gens = [S, T, Tinv]
    h = S.shape[1]
    v = np.broadcast_to(v.astype(complex), (S.shape[0], h)).copy()
    key = lambda u: np.round(np.concatenate([u.real.ravel(), u.imag.ravel()]) * 1e5).astype(np.int64).tobytes()
    seen = {key(v): 0}; vecs = [v]; frontier = [v]
    while frontier and len(vecs) < maxn:
        new = []
        for u in frontier:
            for g in gens:
                w = np.einsum('aij,aj->ai', g, u)
                k = key(w)
                if k not in seen:
                    seen[k] = len(vecs); vecs.append(w); new.append(w)
        frontier = new
    return np.array(vecs)

def realify(vecs, p):
    """rows: real vectors with <.,.> = B; also masks for QR embeddings."""
    n, h, _ = vecs.shape
    X = np.concatenate([vecs.real, vecs.imag], axis=2) * math.sqrt(2)   # (n, h, 2h)
    qr = np.array([legendre(a, p) == 1 for a in range(1, h + 1)])
    Xp = X[:, qr, :].reshape(n, -1); Xm = X[:, ~qr, :].reshape(n, -1)
    return Xp, Xm

def hnf_rows(M, mod=None):
    """row-style Hermite reduction of integer matrix (list of lists of Python ints); returns nonzero rows."""
    M = [list(r) for r in M]; nr = len(M); nc = len(M[0]); r0 = 0
    for c in range(nc):
        piv = None
        for r in range(r0, nr):
            if M[r][c] != 0:
                if piv is None: piv = r
                else:
                    a, b = M[piv][c], M[r][c]
                    # extended gcd
                    x0, x1, y0, y1 = 1, 0, 0, 1; aa, bb = a, b
                    while bb != 0:
                        q = aa // bb
                        aa, bb = bb, aa - q * bb
                        x0, x1 = x1, x0 - q * x1; y0, y1 = y1, y0 - q * y1
                    g = aa; u, v = x0, y0
                    rp = [u * M[piv][j] + v * M[r][j] for j in range(nc)]
                    rr = [(-b // g) * M[piv][j] + (a // g) * M[r][j] for j in range(nc)]
                    if mod is not None:
                        rp = [rp[j] if j == c else rp[j] % mod for j in range(nc)]
                        rr = [x % mod for x in rr]
                    M[piv], M[r] = rp, rr
        if piv is not None:
            M[r0], M[piv] = M[piv], M[r0]; r0 += 1
    return [r for r in M if any(x != 0 for x in r)]

def exact_solve(G, c):
    """solve G x = c exactly; G integer matrix (numpy float, integral), c integer vector. Fractions."""
    n = G.shape[0]
    A = [[Fraction(int(round(G[i, j]))) for j in range(n)] + [Fraction(int(round(c[i])))] for i in range(n)]
    for i in range(n):
        pv = next(r for r in range(i, n) if A[r][i] != 0)
        A[i], A[pv] = A[pv], A[i]
        inv = 1 / A[i][i]
        A[i] = [a * inv for a in A[i]]
        for r in range(n):
            if r != i and A[r][i] != 0:
                f = A[r][i]
                A[r] = [a - f * b for a, b in zip(A[r], A[i])]
    return [A[i][n] for i in range(n)]

def lattice_basis(Xp, Xm, verbose=True, init=None):
    X = np.concatenate([Xp, Xm], axis=1)
    import scipy.linalg as sla
    if init is None:
        Q, R, piv = sla.qr(X.T, pivoting=True, mode='economic')
        d = np.abs(np.diag(R)); rank = int((d > 1e-6 * d[0]).sum())
        Bas = X[list(piv[:rank])].copy()
    else:
        Bas = init.copy(); rank = Bas.shape[0]
    def reduce(Bas):
        G = np.round(Bas @ Bas.T)
        G2, U = V.lll_gram(G)
        return (U.T @ Bas)
    Bas = reduce(Bas)
    rounds = 0
    while rounds < 400:
        rounds += 1
        G = np.round(Bas @ Bas.T)
        C = np.linalg.solve(G, Bas @ X.T).T
        frac = np.abs(C - np.round(C)).max(axis=1)
        bad = np.where(frac > 1e-6)[0]
        if len(bad) == 0: break
        q = None; x = None
        for j in bad[:300]:
            xx = C[j]; qq = 1; ok = True
            for t in xx:
                f = Fraction(float(t)).limit_denominator(100000)
                if abs(float(f) - t) > 1e-9: ok = False; break
                qq = qq * f.denominator // math.gcd(qq, f.denominator)
                if qq > 10**7: ok = False; break
            if ok and np.abs(qq * xx - np.round(qq * xx)).max() < 1e-5 and (q is None or qq < q):
                q, x = qq, [int(round(qq * t)) for t in xx]
                if q <= 3: break
        if q is None:
            # exact fallback on the first bad generator
            j = bad[0]
            xe = exact_solve(G, np.round(Bas @ X[j]))
            q = 1
            for f in xe: q = q * f.denominator // math.gcd(q, f.denominator)
            x = [int(f * q) for f in xe]
            if verbose: print('   (exact solve used; q has %d digits)' % len(str(q)), flush=True)
        rows = [[q if i == k else 0 for k in range(rank)] for i in range(rank)] + [[a % q for a in x]]
        Hn = hnf_rows(rows, mod=q)
        Hn = [r for r in Hn if any(v != 0 for v in r)]
        # complete: the lattice q Z^n + Z a ; make sure we have rank rows (add q e_j rows if modular reduction dropped some)
        if len(Hn) != rank:
            Hn = hnf_rows(Hn + [[q if i == k else 0 for k in range(rank)] for i in range(rank)])
        assert len(Hn) == rank, (len(Hn), rank)
        Hf = np.array([[float(Fraction(v, q)) for v in r] for r in Hn])
        Bas = Hf @ Bas
        Bas = reduce(Bas)
        if verbose:
            sg, ld = np.linalg.slogdet(np.round(Bas @ Bas.T))
            print('   enlarged (index divides q=%s) logdet now %.4f' % (str(q) if q < 10**12 else '~1e%d' % (len(str(q)) - 1), ld), flush=True)
    return Bas, rank

def run(p, do_enum=True):
    t0 = time.time()
    h = (p - 1) // 2; n = 2 * (p - 1)
    for kappa in (1, -1):
        S, T = build(p, kappa)
        rel = check_relations(S, T, p)
        if rel['S^2=-I'] < 1e-9 and rel['(ST)^3=-I'] < 1e-9: break
    print('p', p, 'kappa', kappa, rel)
    v = v0_vec(p)
    O = orbit(S, T, v)
    print(' orbit size', len(O), 'expected <=', p * (p * p - 1), 'time', round(time.time() - t0, 1))
    Xp, Xm = realify(O, p)
    Bas, rank = lattice_basis(Xp, Xm)
    print(' rank', rank, 'expected', n)
    dp = Xp.shape[1]
    Gp = Bas[:, :dp] @ Bas[:, :dp].T; Gm = Bas[:, dp:] @ Bas[:, dp:].T
    B = Gp + Gm; E = math.sqrt(p) * (Gp - Gm)
    print(' integrality of B, E:', float(np.abs(B - np.round(B)).max()), float(np.abs(E - np.round(E)).max()))
    Bi = np.round(B); Ei = np.round(E)
    print(' E B^-1 E = p B residual', float(np.abs(Ei @ np.linalg.solve(Bi, Ei) - p * Bi).max()))
    sign, logdet = np.linalg.slogdet(Bi)
    D = math.exp(logdet / (p - 1)); theta = (p - 1) ** 2 / (4 * p)
    print(' det B = exp(%.6f) ~ %.6g ; D = det^(1/(p-1)) = %.6f ; theta_p = %.6f' % (logdet, math.exp(logdet), D, theta))
    np.save('B_p%d.npy' % p, Bi); np.save('E_p%d.npy' % p, Ei)
    res = {'p': p, 'rank': rank, 'logdet': logdet, 'D': D, 'theta': theta}
    if do_enum:
        eps = {5: (1 + 5 ** .5) / 2, 13: (3 + 13 ** .5) / 2, 17: 4 + 17 ** .5, 29: (5 + 29 ** .5) / 2, 37: 6 + 37 ** .5, 61: (39 + 5 * 61 ** .5) / 2}[p]
        L = math.log(eps); best = (1e99, None, None)
        ts = np.arange(-L, L + 1e-9, 0.1)
        bound0 = 2 * math.sqrt(theta * D) * math.cosh(0.05) * (1 + 1e-9)
        for t in ts:
            Gt = math.exp(t) * Gp + math.exp(-t) * Gm
            try:
                Xs, vals, flag = V.short_vectors(Gt, bound0 * 1.0000001, maxc=200000, maxnodes=200_000_000)
            except RuntimeError as e:
                print('  t=%.2f enumeration limit' % t); continue
            if len(Xs):
                mp = np.einsum('ij,jk,ik->i', Xs, Gp, Xs); mm = np.einsum('ij,jk,ik->i', Xs, Gm, Xs)
                gam = mp * mm / D
                i = int(np.argmin(gam))
                if gam[i] < best[0]: best = (float(gam[i]), float(t), Xs[i].copy())
                print('  t=%+.2f  #vectors with Q_t<=bound: %d   min Gamma %.6f  (theta %.6f)' % (t, len(Xs), gam.min(), theta))
            else:
                print('  t=%+.2f  none below bound' % t)
        print(' RESULT p=%d: Gamma_min found = %s vs theta = %.6f ; ratio %.6f' % (p, best[0], theta, best[0] / theta if best[1] is not None else float('nan')))
        res['gamma_min'] = best[0]
    print(' total time', round(time.time() - t0, 1))
    return res

if __name__ == '__main__':
    for p in [int(a) for a in sys.argv[1:]]:
        run(p)
