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
        T[ia] = np.diag(np.exp(2j * np.pi * ((a * xs * xs) % p) / p))
        e = (2 * a * np.outer(xs, xs)) % p
        S[ia] = kappa * (2j * np.sin(2 * np.pi * e / p)) / (legendre(a, p) * math.sqrt(p))
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

def lattice_basis_exact(Y, verbose=False):
    """Z-basis of the lattice generated by the rows of Y (real vectors whose Gram matrix is integral),
    using exact rational arithmetic on the integer Gram matrix. Returns (Bas, rank)."""
    import scipy.linalg as sla
    K = Y.shape[0]
    Gall = Y @ Y.T
    Gi = np.round(Gall)
    assert np.abs(Gall - Gi).max() < 1e-3, 'Gram not integral: %g' % np.abs(Gall - Gi).max()
    # choose a well conditioned independent subset: pivoted Cholesky on the integer Gram
    Q, Rm, piv = sla.qr(Y.T, pivoting=True, mode='economic')
    d = np.abs(np.diag(Rm)); rank = int((d > 1e-6 * d[0]).sum())
    idx = [int(i) for i in piv[:rank]]
    Gb = [[Fraction(int(Gi[i, j])) for j in idx] for i in idx]
    # exact inverse times G_{b,all}
    n = rank
    A = [Gb[i] + [Fraction(int(Gi[idx[i], j])) for j in range(K)] for i in range(n)]
    for i in range(n):
        pv = next(rr for rr in range(i, n) if A[rr][i] != 0)
        A[i], A[pv] = A[pv], A[i]
        inv = 1 / A[i][i]
        A[i] = [a * inv for a in A[i]]
        for rr in range(n):
            if rr != i and A[rr][i] != 0:
                f = A[rr][i]
                A[rr] = [a - f * b for a, b in zip(A[rr], A[i])]
    X = [[A[i][n + j] for i in range(n)] for j in range(K)]       # coordinates of generator j in basis b (exact)
    q = 1
    for row in X:
        for f in row: q = q * f.denominator // math.gcd(q, f.denominator)
    rows = [[int(f * q) % q for f in row] for row in X]
    rows = [r for r in rows if any(r)]
    qI = [[q if i == k else 0 for k in range(n)] for i in range(n)]
    H1 = hnf_rows(qI + rows, mod=q)
    H2 = hnf_rows([r for r in H1 if any(r)] + qI)
    H2 = [r for r in H2 if any(r)]
    assert len(H2) == n, (len(H2), n)
    # full Hermite reduction: make the entries above each pivot small
    pivcol = []
    for i, r in enumerate(H2):
        c = next(j for j in range(n) if r[j] != 0); pivcol.append(c)
        if r[c] < 0: H2[i] = [-x for x in r]
    for i in range(n):
        c = pivcol[i]; d = H2[i][c]
        for k in range(n):
            if k != i and H2[k][c] != 0:
                f = H2[k][c] // d
                if f: H2[k] = [a - f * b for a, b in zip(H2[k], H2[i])]
    # exact Gram of the HNF basis:  H Gb H^T / q^2  (integers)
    Gb_int = [[int(Gi[i, j]) for j in idx] for i in idx]
    HG = [[sum(H2[i][k] * Gb_int[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    Gex = [[sum(HG[i][k] * H2[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
    assert all(v % (q * q) == 0 for row in Gex for v in row), 'new Gram not integral'
    Gz = [[v // (q * q) for v in row] for row in Gex]
    # exact pairwise (greedy) reduction on the integer Gram to remove the HNF skew
    Uz = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    changed = True; it = 0
    while changed and it < 100000:
        changed = False; it += 1
        for i in range(n):
            for j in range(n):
                if i == j: continue
                num, den = Gz[i][j], Gz[j][j]
                k = (2 * num + den) // (2 * den)          # nearest integer to num/den
                if k != 0 and (Gz[i][i] - 2 * k * num + k * k * den) < Gz[i][i]:
                    for c in range(n): Gz[i][c] -= k * Gz[j][c]
                    for c in range(n): Gz[c][i] -= k * Gz[c][j]
                    Uz[i] = [x - k * y for x, y in zip(Uz[i], Uz[j])]
                    changed = True
    # final rational transformation  C = Uz * H2 / q  (rows: new basis in terms of the chosen generators)
    Cz = [[Fraction(sum(Uz[i][k] * H2[k][j] for k in range(n)), q) for j in range(n)] for i in range(n)]
    Cf = np.array([[float(f) for f in row] for row in Cz])
    Gred, Urows = lll_gram_fraction(Gz)
    # total rational transformation  (Urows Uz H2)/q  applied to the chosen generators
    UzT = [[sum(Urows[i][k] * Uz[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    Cz = [[Fraction(sum(UzT[i][k] * H2[k][j] for k in range(n)), q) for j in range(n)] for i in range(n)]
    Cf = np.array([[float(f) for f in row] for row in Cz])
    Bas = Cf @ Y[idx]
    if verbose: print('   exact lattice basis: rank %d, common denominator has %d digits' % (n, len(str(q))))
    return Bas, rank

def lll_gram_int(Gz, delta=0.99, maxit=200000):
    """LLL on an exact integer Gram matrix (list of lists of Python ints). Exact integer updates of G and U;
    floating Cholesky (recomputed from the exact G) only steers the algorithm. Returns (G float, U int ndarray), G = U^T G0 U."""
    n = len(Gz)
    G = np.array(Gz, dtype=object)
    U = np.array([[1 if i == j else 0 for j in range(n)] for i in range(n)], dtype=object)
    def chol():
        return np.linalg.cholesky(np.array(G, dtype=float))
    L = chol(); k = 1; it = 0
    while k < n and it < maxit:
        it += 1
        for j in range(k - 1, -1, -1):
            r = int(round(L[k, j] / L[j, j]))
            if r != 0:
                U[:, k] -= r * U[:, j]
                G[k, :] -= r * G[j, :]
                G[:, k] -= r * G[:, j]
                L[k, :] -= r * L[j, :]
        mu = L[k, k - 1] / L[k - 1, k - 1]
        if L[k, k] ** 2 >= (delta - mu * mu) * L[k - 1, k - 1] ** 2:
            k += 1
        else:
            U[:, [k - 1, k]] = U[:, [k, k - 1]]
            G[[k - 1, k], :] = G[[k, k - 1], :]
            G[:, [k - 1, k]] = G[:, [k, k - 1]]
            L = chol()
            k = max(k - 1, 1)
    return np.array(G, dtype=float), U

def lll_gram_fraction(Gz, delta=Fraction(3, 4)):
    """Textbook LLL (Cohen, Alg. 2.6.3) on an exact integer Gram matrix with rational Gram-Schmidt data.
    Returns (G as list of lists of ints, U as list of lists of ints) with G = U^T G0 U (columns of U = new basis)."""
    n = len(Gz)
    G = [[int(x) for x in row] for row in Gz]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]      # U[i] = coordinates of new basis vector i (rows)
    mu = [[Fraction(0)] * n for _ in range(n)]; Bs = [Fraction(0)] * n
    def gso_row(i):
        for j in range(i):
            s = Fraction(G[i][j])
            for l in range(j): s -= mu[j][l] * mu[i][l] * Bs[l]
            mu[i][j] = s / Bs[j]
        s = Fraction(G[i][i])
        for l in range(i): s -= mu[i][l] * mu[i][l] * Bs[l]
        Bs[i] = s
    for i in range(n): gso_row(i)
    def red(k, j):
        r = round(mu[k][j])
        if r == 0: return
        for c in range(n): G[k][c] -= r * G[j][c]
        for c in range(n): G[c][k] -= r * G[c][j]
        U[k] = [a - r * b for a, b in zip(U[k], U[j])]
        for l in range(j): mu[k][l] -= r * mu[j][l]
        mu[k][j] -= r
    k = 1
    while k < n:
        red(k, k - 1)
        if Bs[k] < (delta - mu[k][k - 1] ** 2) * Bs[k - 1]:
            m = mu[k][k - 1]; B = Bs[k] + m * m * Bs[k - 1]
            mu[k][k - 1] = m * Bs[k - 1] / B
            Bs[k] = Bs[k - 1] * Bs[k] / B; Bs[k - 1] = B
            G[k - 1], G[k] = G[k], G[k - 1]
            for row in G: row[k - 1], row[k] = row[k], row[k - 1]
            U[k - 1], U[k] = U[k], U[k - 1]
            for j in range(k - 1): mu[k - 1][j], mu[k][j] = mu[k][j], mu[k - 1][j]
            for i in range(k + 1, n):
                t = mu[i][k]
                mu[i][k] = mu[i][k - 1] - m * t
                mu[i][k - 1] = t + mu[k][k - 1] * mu[i][k]
            k = max(k - 1, 1)
        else:
            for j in range(k - 2, -1, -1): red(k, j)
            k += 1
    return G, U
