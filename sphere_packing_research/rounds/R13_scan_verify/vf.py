"""Adversarial verifier (R13): independent recomputation of fixed sublattices of the SL_2(p) odd Weil lattice.
Own code: Weil rep, sampling, norm map, exact lattice-from-Gram (Fractions + modular HNF), exact unimodular Gram
Gu = (b p B -/+ a E)/(2 p r m), exact s = (b E[x] - a B[x])/(4 r m), exact Gamma/theta = (p B[x]^2 - E[x]^2)/(r^2 p (p-1)^2)."""
import numpy as np, math, sys, time, collections
from fractions import Fraction
sys.path.insert(0, '/home/user/Mathathon/sphere_packing_research/rounds/R2_voronoi_descent')
import voronoi as V
def leg(a, p):
    a = int(a) % p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
def unit(p):
    b = 1
    while True:
        a = math.isqrt(p * b * b - 4)
        if a * a == p * b * b - 4: return a, b
        b += 1
def sample(p, K, seed):
    h = (p - 1) // 2; r = (p - 1) // 4
    xs = np.arange(1, h + 1); A = np.arange(1, h + 1)
    lg = np.array([leg(a, p) for a in A]); z = np.exp(2j * np.pi * np.arange(p) / p)
    e = (2 * A[:, None, None] * xs[None, :, None] * xs[None, None, :]) % p
    S = -(z[e] - z[(-e) % p]) / (lg[:, None, None] * math.sqrt(p)); del e
    Tex = (A[:, None] * (xs * xs)[None, :]) % p
    I = np.eye(h); T = z[Tex]; ST = S[:3] * T[:3, None, :]
    chk = (float(np.abs(S[:3] @ S[:3] + I).max()), float(np.abs(ST @ ST @ ST + I).max()))
    v0 = np.zeros(h)
    for x in range(1, p):
        if pow(x, r, p) == 1:
            if x <= h: v0[x - 1] += 1
            else: v0[p - x - 1] -= 1
    rng = np.random.default_rng(seed)
    U = np.broadcast_to(v0.astype(complex)[None, :, None], (h, h, K)).copy()
    wl = rng.integers(5, 10, size=K)
    for step in range(9):
        k = rng.integers(0, p, size=K)
        Un = np.matmul(S, U * z[(Tex[:, :, None] * k[None, None, :]) % p])
        U = np.where((step < wl)[None, None, :], Un, U)
    k = rng.integers(0, p, size=K)
    U = U * z[(Tex[:, :, None] * k[None, None, :]) % p]
    return U, lg, chk
def normmat(p, m):
    h = (p - 1) // 2
    fac = [q for q in range(2, p) if (p - 1) % q == 0 and all(q % d for d in range(2, q))]
    g0 = next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in fac))
    c = pow(g0, (p - 1) // m, p); N = np.zeros((h, h))
    for k in range(m):
        ck = pow(c, k, p); s = leg(ck, p)
        for x in range(1, h + 1):
            y = (ck * x) % p
            if y <= h: N[x - 1, y - 1] += s
            else: N[x - 1, p - y - 1] -= s
    return N
def xgcd(a, b):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q = a // b; a, b = b, a - q * b; x0, x1 = x1, x0 - q * x1; y0, y1 = y1, y0 - q * y1
    return a, x0, y0
def detint(M):
    n = len(M); A = [[Fraction(int(x)) for x in r] for r in M]; d = Fraction(1)
    for i in range(n):
        pv = next((r for r in range(i, n) if A[r][i] != 0), None)
        if pv is None: return 0
        if pv != i: A[i], A[pv] = A[pv], A[i]; d = -d
        d *= A[i][i]
        for r in range(i + 1, n):
            if A[r][i] != 0:
                f = A[r][i] / A[i][i]; A[r] = [x - f * y for x, y in zip(A[r], A[i])]
    return d
def lattice_from_gram(Bg, k):
    """Bg: K x K integer Gram (object array) of generators spanning rank k. Returns (idx, H, q): basis rows = H/q in terms of generators idx."""
    K = Bg.shape[0]; R = np.array(Bg, dtype=float); idx = []
    for _ in range(k):                      # greedy pivoted Cholesky
        d = np.diag(R).copy(); d[idx] = -1; i = int(np.argmax(d)); idx.append(i)
        R = R - np.outer(R[:, i], R[i, :]) / R[i, i]
    resid = float(np.abs(np.diag(R)).max())
    A = [[Fraction(int(Bg[i, j])) for j in idx] + [Fraction(int(Bg[i, j])) for j in range(K)] for i in idx]
    for i in range(k):
        pv = next(rr for rr in range(i, k) if A[rr][i] != 0); A[i], A[pv] = A[pv], A[i]
        inv = 1 / A[i][i]; A[i] = [a * inv for a in A[i]]
        for rr in range(k):
            if rr != i and A[rr][i] != 0:
                f = A[rr][i]; A[rr] = [a - f * b for a, b in zip(A[rr], A[i])]
    q = 1
    for i in range(k):
        for j in range(K): q = q * A[i][k + j].denominator // math.gcd(q, A[i][k + j].denominator)
    rows = [[int(A[i][k + j] * q) % q for i in range(k)] for j in range(K)]
    rows = [r for r in rows if any(r)]; H = []
    for c in range(k):
        piv = [0] * k; piv[c] = q; new = []
        for r in rows:
            if r[c] == 0: new.append(r); continue
            a, b = piv[c], r[c]; g, u, v = xgcd(a, b)
            npv = [u * piv[j] + v * r[j] for j in range(k)]; nr = [((a // g) * r[j] - (b // g) * piv[j]) % q for j in range(k)]
            piv = [npv[j] if j <= c else npv[j] % q for j in range(k)]
            if any(nr): new.append(nr)
        H.append(piv); rows = new
    for c in range(k):                      # reduce above the diagonal
        for i in range(c):
            f = H[i][c] // H[c][c]
            if f: H[i] = [x - f * y for x, y in zip(H[i], H[c])]
    return idx, H, q, resid
def tr(H, G, q):
    k = len(H); HG = [[sum(H[i][l] * int(G[l][j]) for l in range(k)) for j in range(k)] for i in range(k)]
    M = [[sum(HG[i][l] * H[j][l] for l in range(k)) for j in range(k)] for i in range(k)]
    assert all(x % (q * q) == 0 for r in M for x in r), 'non-integral transformed Gram'
    return np.array([[x // (q * q) for x in r] for r in M], dtype=object)
def pair_reduce(BM, EM):
    k = len(BM); B = [[int(x) for x in r] for r in BM]; E = [[int(x) for x in r] for r in EM]; changed = True
    while changed:
        changed = False
        for i in range(k):
            for j in range(k):
                if i == j: continue
                num, den = B[i][j], B[j][j]; c = (2 * num + den) // (2 * den)
                if c != 0 and B[i][i] - 2 * c * num + c * c * den < B[i][i]:
                    for M in (B, E):
                        M[i] = [x - c * y for x, y in zip(M[i], M[j])]
                        for rr in range(k): M[rr][i] -= c * M[rr][j]
                    changed = True
    return np.array(B, dtype=object), np.array(E, dtype=object)
def analyse(p, m, Bg, Eg, tag='', do4=True, tscan=True):
    t0 = time.time(); r = (p - 1) // 4; k = 2 * (p - 1) // m; a, b = unit(p); n = 2 * (p - 1)
    idx, H, q, resid = lattice_from_gram(Bg, k)
    Bi = [[int(Bg[i, j]) for j in idx] for i in idx]; Ei = [[int(Eg[i, j]) for j in idx] for i in idx]
    BM = tr(H, Bi, q); EM = tr(H, Ei, q); BM, EM = pair_reduce(BM, EM)
    sys.path.insert(0, '/home/user/Mathathon/sphere_packing_research/rounds/R10_weil_lattices/orchestrator'); import proto as P
    Gred, Ur = P.lll_gram_fraction([[int(x) for x in r] for r in BM]); Uo = np.array(Ur, dtype=object).T   # exact LLL (reused routine; result re-verified exactly below)
    BM = Uo.T.dot(BM).dot(Uo); EM = Uo.T.dot(EM).dot(Uo)
    detB = detint(BM); expect = Fraction(m) ** k * Fraction(r * r * p) ** (k // 2)
    print('%s p=%d C_%d rank %d: rank-resid %.1e, q digits %d, det_B(M)/(m^k (r^2 p)^(k/2)) = %s  [exact; =1 means log R_U=0]' % (tag, p, m, k, resid, len(str(q)), detB / expect), flush=True)
    # O_K stability: (B+E)/2 B^-1 integral ?
    Bf = np.array(BM, dtype=float); half = (BM + EM)
    par = all(int(x) % 2 == 0 for x in half.ravel())
    C = np.linalg.solve(Bf.T, (np.array(half, dtype=float) / 2).T).T; Ci = np.round(C).astype(np.int64).astype(object)
    okst = par and bool((Ci.dot(BM) * 2 == half).all())
    C2 = np.linalg.solve(Bf.T, np.array(EM, dtype=float).T).T; C2i = np.round(C2).astype(np.int64).astype(object)
    sqst = bool((C2i.dot(BM) == EM).all())
    print('   Z[sqrt p]-stable: %s ; O_K-stable: %s (max frac part of omega-matrix %.2e)' % (sqst, okst, float(np.abs(C - np.round(C)).max())), flush=True)
    out = {'detratio': detB / expect, 'okst': okst}
    for sgn in (-1, 1):
        num = b * p * BM + sgn * a * EM; den = 2 * p * r * m
        if all(int(x) % den == 0 for x in num.ravel()):
            Gu = np.array([[int(x) // den for x in row] for row in num], dtype=object)
            du = detint(Gu); even = all(int(Gu[i, i]) % 2 == 0 for i in range(k))
            print('   unimodular member (sign %+d): EXACTLY integral, even %s, det = %s' % (sgn, even, du), flush=True)
            Guf = np.array(Gu, dtype=float)
            Xs, vals, flag = V.short_vectors(Guf, 4.01 if do4 else 2.01, maxc=400000, maxnodes=2_000_000_000)
            X = Xs.astype(object)
            Nn = np.array([int(x.dot(Gu).dot(x)) for x in X]); Bx = np.array([int(x.dot(BM).dot(x)) for x in X], dtype=object); Ex = np.array([int(x.dot(EM).dot(x)) for x in X], dtype=object)
            snum = b * Ex + sgn * a * Bx * (1 if sgn == -1 else -1) if False else (b * Ex - a * Bx if sgn == -1 else b * Ex + a * Bx)
            sint = all(int(x) % (4 * r * m) == 0 for x in snum); s = [abs(Fraction(int(x), 4 * r * m)) for x in snum]
            cnt = collections.Counter(zip(Nn.tolist(), [str(x) for x in s]))
            gam = [Fraction(int(p * bx * bx - ex * ex), r * r * p * (p - 1) ** 2) for bx, ex in zip(Bx, Ex)]
            print('   enum flag %d; #roots %d, #norm4 %d ; s integral %s ; min Gamma/theta over norm<=4: %s' % (flag, 2 * int((Nn == 2).sum()), 2 * int((Nn == 4).sum()), sint, min(gam) if gam else None), flush=True)
            print('   (norm,|s|): pairs ->', sorted(cnt.items(), key=lambda t: (t[0][0], float(Fraction(t[0][1])))), flush=True)
            Rt = Xs[Nn == 2]
            if len(Rt):
                Rts = np.concatenate([Rt, -Rt]); IP = np.round(Rts @ Guf @ Rts.T).astype(int); nR = len(Rts); seen = np.zeros(nR, bool); comps = []
                for i in range(nR):
                    if seen[i]: continue
                    st = [i]; seen[i] = True; c = [i]
                    while st:
                        x = st.pop()
                        for y in np.where((IP[x] != 0) & (~seen))[0]: seen[y] = True; st.append(y); c.append(y)
                    comps.append(c)
                print('   root components (roots, rank):', sorted((len(c), int(np.linalg.matrix_rank(IP[np.ix_(c, c)].astype(float)))) for c in comps), flush=True)
            out.update({'detGu': du, 'roots': 2 * int((Nn == 2).sum()), 'n4': 2 * int((Nn == 4).sum())})
            break
    else:
        print('   NO exactly integral unimodular-normalised member at t=+-log eps', flush=True)
    if tscan:
        L = math.log((a + b * math.sqrt(p)) / 2); D = r * r * p; theta = (p - 1) ** 2 / (4 * p); step = 0.1
        Gp = (Bf + np.array(EM, dtype=float) / math.sqrt(p)) / 2; Gm = (Bf - np.array(EM, dtype=float) / math.sqrt(p)) / 2
        bound = 2 * math.sqrt(theta * D) * math.cosh(step / 2) * (1 + 1e-9); best = None; tot = 0; fail = 0; f2 = 0
        for t in np.arange(-L, L + step, step):
            try: Xs, vals, flag = V.short_vectors(math.exp(t) * Gp + math.exp(-t) * Gm, bound, maxc=200000, maxnodes=500_000_000)
            except Exception as ex: fail += 1; continue
            f2 += (flag == 2); tot += len(Xs)
            for x in Xs.astype(object):
                bx = int(x.dot(BM).dot(x)); ex = int(x.dot(EM).dot(x)); g = Fraction(p * bx * bx - ex * ex, r * r * p * (p - 1) ** 2)
                if best is None or g < best: best = g
        print('   t-circle scan (step %.2f, %d grid pts incl. wrap): hits %d, failures %d, cap-hits %d, min Gamma/theta = %s (exact)  [%.1fs]' % (step, len(np.arange(-L, L + step, step)), tot, fail, f2, best, time.time() - t0), flush=True)
    return out
if __name__ == '__main__':
    p = int(sys.argv[1]); ms = [int(x) for x in sys.argv[2].split(',')]; K = int(sys.argv[3]); seed = int(sys.argv[4]); do4 = (len(sys.argv) <= 5 or sys.argv[5] != 'no4')
    t0 = time.time(); U, lg, chk = sample(p, K, seed); h = (p - 1) // 2
    print('p=%d K=%d seed=%d: relation residuals (S^2+I, (ST)^3+I on 3 embeddings) %s, sampling %.1fs' % (p, K, seed, chk, time.time() - t0), flush=True)
    for m in ms:
        N = normmat(p, m); W = np.matmul(N, U)
        Wp = W[lg == 1].reshape(-1, K); Wm = W[lg == -1].reshape(-1, K)
        Gp = 2 * (Wp.conj().T @ Wp).real; Gm = 2 * (Wm.conj().T @ Wm).real
        Bf = Gp + Gm; Ef = math.sqrt(p) * (Gp - Gm)
        print(' C_%d: float integrality of B %.1e, of E %.1e' % (m, np.abs(Bf - np.round(Bf)).max(), np.abs(Ef - np.round(Ef)).max()), flush=True)
        Bg = np.round(Bf).astype(np.int64).astype(object); Eg = np.round(Ef).astype(np.int64).astype(object)
        for K0 in ([40, 160, K] if len(sys.argv) > 6 else [K]):
            res = analyse(p, m, Bg[:K0, :K0], Eg[:K0, :K0], tag='[first %d samples]' % K0, do4=do4, tscan=(K0 == K))
        if not res['okst'] or len(sys.argv) > 6:
            # O_K-closure: generators g_i and omega g_i, omega=(1+sqrt p)/2
            Hh = (Bg + Eg) // 2; c = (p - 1) // 4
            Bx = np.block([[Bg, Hh], [Hh, Hh + c * Bg]]); PE = (Eg + p * Bg) // 2; Ex = np.block([[Eg, PE], [PE, PE + c * Eg]])
            analyse(p, m, Bx, Ex, tag='[O_K-closure of all %d samples]' % K, do4=do4, tscan=False)
    print('total %.1fs' % (time.time() - t0), flush=True)
