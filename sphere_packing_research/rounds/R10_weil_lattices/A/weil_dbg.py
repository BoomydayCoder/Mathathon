#!/usr/bin/env python3
"""Weil-representation lattices for SL_2(p), p = 1 mod 4: build Lambda_0 = Z[G] v0, order check, Gamma test."""
import numpy as np, math, sys, os, json, time, ctypes
import scipy.linalg
HERE = os.path.dirname(os.path.abspath(__file__))
lat = ctypes.CDLL(os.path.join(HERE, 'lat.so'))
fpl = ctypes.CDLL(os.path.join(HERE, 'fp_enum.so'))
fpl.fp_enum.restype = ctypes.c_long
T0 = time.time()
def log(*a):
    print(f'[{time.time()-T0:7.1f}s]', *a, flush=True)

def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)

def fund_unit(p):
    b = 1
    while True:
        for s in (-4, 4):
            a2 = p * b * b + s
            a = math.isqrt(a2)
            if a * a == a2:
                return a, b, (a + b * math.sqrt(p)) / 2
        b += 1

def xgcd(a, b):
    s0, s1, t0, t1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    if a < 0:
        a, s0, t0 = -a, -s0, -t0
    return a, s0, t0

def build(p, kappa=1):
    h = (p - 1) // 2
    x = np.arange(1, h + 1)
    xy = np.outer(x, x)
    S = np.zeros((h, h, h), complex)
    T = np.zeros((h, h, h), complex)
    sq = math.sqrt(p)
    for ai in range(h):
        a = ai + 1
        e = lambda k: np.exp(2j * np.pi * ((a * k) % p) / p)
        S[ai] = kappa * (e(2 * xy) - e(-2 * xy)) / (legendre(a, p) * sq)
        T[ai] = np.diag(e(x * x))
    return S, T

def mm(A, B):
    return np.einsum('axy,ayz->axz', A, B)

def act(g, U):
    return np.einsum('axy,...ay->...ax', g, U)

def torus(p, a):
    h = (p - 1) // 2
    Dm = np.zeros((h, h))
    for x in range(1, h + 1):
        y = (a * x) % p
        if y <= h:
            Dm[y - 1, x - 1] = 1
        else:
            Dm[p - y - 1, x - 1] = -1
    return np.tile(Dm.astype(complex), (h, 1, 1))

class Ctx:
    pass

def remb(U, c):
    sh = U.shape[:-2]
    P = U[..., c.plus, :].reshape(sh + (-1,))
    M = U[..., c.minus, :].reshape(sh + (-1,))
    return math.sqrt(2) * np.concatenate([P.real, P.imag, M.real, M.imag], axis=-1)

def orbit(V0, gens, limit=10**7):
    def key(V):
        return np.round(np.concatenate([V[0].real, V[0].imag]) * 1e5).astype(np.int64).tobytes()
    vecs = [V0]
    seen = {key(V0): 0}
    par = [(-1, -1)]
    frontier = [0]
    while frontier:
        F = np.array([vecs[i] for i in frontier])
        new = []
        for gi, g in enumerate(gens):
            GF = np.einsum('axy,nay->nax', g, F)
            for k, i in enumerate(frontier):
                kk = key(GF[k])
                if kk not in seen:
                    seen[kk] = len(vecs)
                    vecs.append(GF[k])
                    par.append((i, gi))
                    new.append(len(vecs) - 1)
        frontier = new
        if len(vecs) > limit:
            break
    return np.array(vecs), par, seen, key

def c_lll(G, delta=0.99):
    n = len(G)
    L = np.ascontiguousarray(np.linalg.cholesky((G + G.T) / 2))
    U = np.ascontiguousarray(np.eye(n, dtype=np.int64))
    mu = np.zeros((n, n)); r = np.zeros(n)
    rc = lat.lll_c(n, n, L.ctypes.data, U.ctypes.data, ctypes.c_double(delta), mu.ctypes.data, r.ctypes.data)
    assert rc == 0
    return U

def c_bkz(G, beta, tours=8, maxnodes=20_000_000, U0=None, delta=0.99):
    n = len(G)
    L = np.ascontiguousarray(np.linalg.cholesky((G + G.T) / 2))
    U = np.ascontiguousarray(np.eye(n, dtype=np.int64) if U0 is None else U0.copy())
    r = np.zeros(n)
    nt = lat.bkz_c(n, n, L.ctypes.data, U.ctypes.data, int(beta), int(tours), ctypes.c_double(delta),
                   ctypes.c_long(maxnodes), r.ctypes.data)
    return U, r, nt

def enum_short(G, bound, maxc=400000, maxnodes=4_000_000_000):
    """all nonzero x (one per +-pair) with x G x^T <= bound; rows x in the given basis."""
    n = len(G)
    U = c_lll(G)
    G2 = U @ G @ U.T
    L = np.ascontiguousarray(np.linalg.cholesky((G2 + G2.T) / 2))
    out = np.zeros((maxc, n), dtype=np.int32); vals = np.zeros(maxc)
    flag = ctypes.c_int(0); nodes = ctypes.c_long(0)
    cnt = fpl.fp_enum(n, L.ctypes.data, ctypes.c_double(bound), ctypes.c_long(maxc), out.ctypes.data, vals.ctypes.data,
                      ctypes.c_long(maxnodes), ctypes.byref(flag), ctypes.byref(nodes))
    Y = out[:cnt].astype(np.int64)
    return Y @ U, flag.value, nodes.value

def add_vec_basis(q, f):
    """triangular basis (rows) of q Z^n + Z f."""
    n = len(f)
    P = np.zeros((n, n), dtype=object)
    for i in range(n):
        P[i, i] = q
    v = np.array([int(z) % q for z in f], dtype=object)
    for j in range(n):
        vj = v[j]
        if vj == 0:
            continue
        g, s, t = xgcd(q, vj)
        newp = s * P[j] + t * v
        v = (q // g) * v - (vj // g) * P[j]
        newp[j + 1:] = newp[j + 1:] % q
        v[j + 1:] = v[j + 1:] % q
        P[j] = newp
    assert all(z == 0 for z in v)
    for j in range(n - 1, -1, -1):
        for k in range(j + 1, n):
            cq = P[j, k] // P[k, k]
            if cq:
                P[j] = P[j] - cq * P[k]
    return P

def bareiss_det(M):
    A = np.array([[int(z) for z in row] for row in M], dtype=object)
    n = len(A); sign = 1; prev = 1
    for k in range(n - 1):
        if A[k, k] == 0:
            nz = [i for i in range(k + 1, n) if A[i, k] != 0]
            if not nz:
                return 0
            A[[k, nz[0]]] = A[[nz[0], k]]; sign = -sign
        A[k + 1:, k + 1:] = (A[k + 1:, k + 1:] * A[k, k] - np.outer(A[k + 1:, k], A[k, k + 1:])) // prev
        prev = A[k, k]
    return sign * A[n - 1, n - 1]

def saturate(Y):
    Yo = np.array([[int(z) for z in row] for row in Y], dtype=object)
    r, n = Yo.shape
    A = Yo.copy()
    for i in range(r):
        for j in range(i + 1, n):
            if A[i, j] != 0:
                a, b = A[i, i], A[i, j]
                g, s, t = xgcd(a, b)
                ci = A[:, i].copy(); cj = A[:, j].copy()
                A[:, i] = s * ci + t * cj
                A[:, j] = -(b // g) * ci + (a // g) * cj
        assert A[i, i] != 0
    Tm = A[:, :r]
    F = np.zeros((r, n), dtype=object)
    for i in range(r):
        acc = Yo[i].copy()
        for k in range(i):
            acc = acc - Tm[i, k] * F[k]
        d = Tm[i, i]
        assert all(z % d == 0 for z in acc)
        F[i] = acc // d
    return F.astype(np.int64), Tm

def setup(p, full_orbit=True):
    c = Ctx(); c.p = p; h = c.h = (p - 1) // 2; c.n = 4 * h
    S, T = build(p, 1)
    I = np.tile(np.eye(h, dtype=complex), (h, 1, 1))
    c.I = I
    e2 = np.abs(mm(S, S) + I).max()
    ST = mm(S, T); ST3 = mm(ST, mm(ST, ST))
    r3p = np.abs(ST3 - I).max(); r3m = np.abs(ST3 + I).max()
    Tp = I.copy()
    for _ in range(p):
        Tp = mm(Tp, T)
    unit = max(np.abs(mm(S, S.conj().transpose(0, 2, 1)) - I).max(), np.abs(mm(T, T.conj().transpose(0, 2, 1)) - I).max())
    log(f'p={p} relations: |S^2+I|={e2:.1e} |(ST)^3-I|={r3p:.1e} |(ST)^3+I|={r3m:.1e} |T^p-I|={np.abs(Tp-I).max():.1e} unitarity={unit:.1e}')
    c.S, c.T = S, T
    sg = np.array([legendre(a, p) for a in range(1, h + 1)])
    c.plus = np.where(sg == 1)[0]; c.minus = np.where(sg == -1)[0]
    c.hp = h * h
    # pure vector
    if p % 8 == 5 and p > 5:
        m = (p - 1) // 4
        mu = [x for x in range(1, p) if pow(x, m, p) == 1]
        v0 = np.zeros(h)
        for x in mu:
            if x <= h: v0[x - 1] += 1
            else: v0[p - x - 1] -= 1
        V0 = np.tile(v0.astype(complex), (h, 1))
        # generator of mu_m
        gen = [x for x in mu if all(pow(x, m // q, p) != 1 for q in range(2, m + 1) if m % q == 0)][0]
        c.Hgen = torus(p, gen); c.Hord = m; c.torus_gen = gen
    elif p == 5:
        V0 = np.zeros((h, h), complex); V0[:, 0] = 1
        c.Hgen = I.copy(); c.Hord = 1
    else:
        g3 = mm(ST, ST)
        assert np.abs(mm(g3, mm(g3, g3)) - I).max() < 1e-9
        e1 = np.zeros((h, h), complex); e1[:, 0] = 1
        V0 = e1 + act(g3, e1) + act(g3, act(g3, e1))
        c.Hgen = g3; c.Hord = 3
    c.V0 = V0
    assert np.abs(act(c.Hgen, V0) - V0).max() < 1e-9
    t1 = time.time()
    OV, par, seen, key = orbit(V0, [S, T])
    c.OV = OV; c.seen = seen; c.key = key
    log(f'orbit size {len(OV)} (|G|={p*(p*p-1)}, |G|/|H|={p*(p*p-1)//c.Hord})')
    # check H generator maps orbit to orbit (sample)
    smp = OV[:: max(1, len(OV) // 200)]
    ok = all(key(v) in seen for v in act(c.Hgen, smp))
    log('H-generator maps orbit sample into orbit:', ok)
    return c

def gauss_jordan_adj(G):
    n = len(G)
    A = np.zeros((n, 2 * n), dtype=object)
    for i in range(n):
        for j in range(n):
            A[i, j] = int(G[i, j])
        A[i, n + i] = 1
    prev = 1
    for k in range(n):
        piv = A[k, k]
        assert piv != 0
        rowk = A[k].copy(); colk = A[:, k].copy()
        A = (piv * A - np.outer(colk, rowk)) // prev
        A[k] = rowk
        prev = piv
    return A[0, 0], A[:, n:]

def hnf_mod(vecs, d, n, P=None):
    """triangular basis of d Z^n + sum Z v (exact integers)."""
    if P is None:
        P = np.zeros((n, n), dtype=object)
        for i in range(n):
            P[i, i] = d
    for v in vecs:
        v = np.array([int(z) % d for z in v], dtype=object)
        for j in range(n):
            vj = v[j]
            if vj == 0:
                continue
            pj = P[j, j]
            if vj % pj == 0:
                v = v - (vj // pj) * P[j]
            else:
                g, s, t = xgcd(pj, vj)
                newp = s * P[j] + t * v
                v = (pj // g) * v - (vj // g) * P[j]
                newp[j + 1:] = newp[j + 1:] % d
                P[j] = newp
            v[j + 1:] = v[j + 1:] % d
    return P

def hnf_reduce(P):
    P = P.copy(); n = len(P)
    for j in range(n - 1, -1, -1):
        for k in range(j + 1, n):
            cq = P[j, k] // P[k, k]
            if cq:
                P[j] = P[j] - cq * P[k]
    return P

def solve_mod(Nmat, rhs, Q):
    """solve y @ Nmat = rhs mod Q (Q prime)."""
    n = len(rhs)
    A = np.empty((n, n + 1), dtype=object)
    A[:, :n] = Nmat.T % Q
    A[:, n] = rhs % Q
    for k in range(n):
        pv = next(i for i in range(k, n) if A[i, k] != 0)
        if pv != k:
            A[[k, pv]] = A[[pv, k]]
        inv = pow(int(A[k, k]), -1, Q)
        rowk = (A[k] * inv) % Q
        colk = A[:, k].copy()
        A = (A - np.outer(colk, rowk)) % Q
        A[k] = rowk
    return A[:, n]

def build_lattice(c, rng):
    p, h, n = c.p, c.h, c.n
    X = remb(c.OV, c)
    N = len(X)
    idx = rng.choice(N, size=min(N, 8 * n), replace=False)
    Q, R, piv = scipy.linalg.qr(X[idx].T, mode='economic', pivoting=True)
    d = np.abs(np.diag(R))
    rank = int((d > 1e-7 * d[0]).sum())
    log(f'rank of orbit sample = {rank} (expected {n}); R diag at rank-1: {d[rank-1]:.3e}, next: {d[rank] if rank < len(d) else 0:.3e}')
    c.rank = rank
    assert rank == n, 'rank mismatch'
    Qn = Q[:, :n]
    Xr = X @ Qn
    spanres = np.abs(X - Xr @ Qn.T).max()
    log(f'all orbit vectors in the {n}-dim span: residual {spanres:.2e}')
    sel = idx[piv[:n]]
    Mr = Xr[sel]
    GMf = Mr @ Mr.T
    GM = np.rint(GMf).astype(np.int64)
    log(f'Gram of M integer residual {np.abs(GMf-GM).max():.2e}')
    dM, adj = gauss_jordan_adj(GM)
    chk = GM.astype(object).dot(adj)
    assert all(chk[i, j] == (dM if i == j else 0) for i in range(n) for j in range(n))
    log(f'det G_M has {len(str(dM))} digits; adjugate verified exactly')
    gen = list(range(min(N, 2 * n))) + list(rng.choice(N, size=min(N, 2 * n), replace=False))
    P = None; rounds = 0
    GMf = GM.astype(float)
    while True:
        IvF = Xr[gen] @ Mr.T
        Iv = np.rint(IvF).astype(np.int64)
        assert np.abs(IvF - Iv).max() < 1e-6
        dc = Iv.astype(object).dot(adj)
        P = hnf_mod(dc, dM, n, P)
        Pr = hnf_reduce(P)
        # chain of prime-index steps from Z^n (M-coordinates) up to the lattice (Pr/dM)
        Ncur = np.zeros((n, n), dtype=object)
        for i in range(n):
            Ncur[i, i] = dM
        Ccur = np.eye(n); steps = []; index = 1
        QQ = 2 ** 127 - 1
        for j in range(n - 1, -1, -1):
            k = dM // Pr[j, j]
            if k == 1:
                continue
            index *= k
            fac = []; kk = k; q = 2
            while kk > 1 and q < 200000:
                while kk % q == 0:
                    fac.append(q); kk //= q
                q += 1 if q == 2 else 2
            if kk > 1:
                fac.append(kk)
            m = k
            for l in fac:
                m //= l
                rhs = np.array([l * m * Pr[j, t] for t in range(n)], dtype=object)
                y = solve_mod(Ncur, rhs, QQ)
                y = np.array([int(z) if int(z) < QQ // 2 else int(z) - QQ for z in y], dtype=object)
                assert all(z == 0 for z in (y.dot(Ncur) - rhs)), 'exact solve failed'
                fi = [int(z) % l for z in y]
                assert any(fi)
                Pl = add_vec_basis(l, fi)
                Cf = np.array([[float(Pl[i, t] / l) if Pl[i, t] else 0.0 for t in range(n)] for i in range(n)]) @ Ccur
                U = c_lll(Cf @ GMf @ Cf.T)
                assert np.abs(U).max() < 2 ** 60
                Nn = U.astype(object).dot(Pl.dot(Ncur))
                assert all(z % l == 0 for z in Nn.ravel())
                Ncur = Nn // l
                Ccur = np.array([[float(Ncur[i, t] / dM) if Ncur[i, t] else 0.0 for t in range(n)] for i in range(n)])
                steps.append(int(l))
                if l > 10**6: log('giant step', l, 'max|U|', np.abs(U).max(), 'detU', bareiss_det(U))
        # exact membership test of the generators in the chain lattice and in HNF lattice
        nbad = 0
        for vv in dc[:10]:
            yv = solve_mod(Ncur, np.array([int(z) for z in vv], dtype=object), QQ)
            yv = np.array([int(z) if int(z) < QQ // 2 else int(z) - QQ for z in yv], dtype=object)
            if any(z != 0 for z in (yv.dot(Ncur) - vv)): nbad += 1
        nbad2 = 0
        for vv in dc[:10]:
            yv = solve_mod(Pr, np.array([int(z) for z in vv], dtype=object), QQ)
            yv = np.array([int(z) if int(z) < QQ // 2 else int(z) - QQ for z in yv], dtype=object)
            if any(z % dM != 0 for z in (yv.dot(Pr) - vv)): nbad2 += 1
        log(f'exact membership: first 10 generators NOT in chain lattice: {nbad}; not in HNF lattice (mod dM): {nbad2}')
        BasC = np.tensordot(Ccur, c.OV[sel], axes=1)
        Br = remb(BasC, c) @ Qn
        C = np.linalg.solve(Br.T, Xr.T).T
        Fr = np.abs(C - np.rint(C)).max(1)
        bad = np.where(Fr > 1e-6)[0]
        rounds += 1
        yv = solve_mod(Ncur, np.array([int(z) for z in dc[0]], dtype=object), QQ)
        yv = np.array([float(int(z) if int(z) < QQ // 2 else int(z) - QQ) for z in yv])
        log('DBG frac of first 10 gens', Fr[gen[:10]], 'exact y vs float C for gen0: max diff', np.abs(yv - C[gen[0]]).max(), 'max|y|', np.abs(yv).max(), 'max|C|', np.abs(C).max())
        cM = np.array([float(int(z) / dM) for z in dc[0]]); cMf = np.linalg.solve(Mr.T, Xr[gen[0]])
        log('DBG M-coords exact vs float: max diff', np.abs(cM - cMf).max(), 'cond(Mr)', np.linalg.cond(Mr), 'cond(Br)', np.linalg.cond(Br), 'cond(Ccur)', np.linalg.cond(Ccur))
        log('DBG |BasC - Ccur@M| check via Br vs Ccur@Mr:', np.abs(Br - Ccur @ Mr).max(), ' solve residual', np.abs(C @ Br - Xr).max())
        GB = Br @ Br.T
        log(f'HNF round {rounds}: generators used {len(gen)}, index so far {index}, prime steps {steps}, orbit vectors not yet integral: {len(bad)}; max|Ccur|={np.abs(Ccur).max():.3e}, Gram diag min/max {GB.diagonal().min():.3f}/{GB.diagonal().max():.3f}, Gram int resid {np.abs(GB-np.rint(GB)).max():.2e}, frac of gens max {Fr[gen].max():.2e}, frac hist {np.histogram(Fr, bins=[0,1e-9,1e-6,1e-3,0.1,0.51])[0]}')
        assert rounds < 6
        if len(bad) == 0:
            break
        gen = list(bad[:80])
    c.index_M = index; c.dM = dM
    Xb = remb(BasC, c)
    U = c_lll(Xb @ Xb.T)
    BasC = np.tensordot(U.astype(float), BasC, axes=1)
    c.BasC = BasC
    finalize(c)
    log(f'EXACT check det(B) * index^2 == det(G_M): {c.detB * c.index_M**2 == c.dM}')
    C = coords(c, c.OV)
    log(f'check: orbit vectors integral in basis, max frac = {np.abs(C-np.rint(C)).max():.2e}')
    for nm, g in (('S', c.S), ('T', c.T)):
        Cg = coords(c, act(g, c.BasC))
        log(f'check: {nm}*basis integral, max frac = {np.abs(Cg-np.rint(Cg)).max():.2e}, |det|={abs(np.linalg.det(np.rint(Cg))):.6f}')

def finalize(c):
    p = c.p
    Xb = remb(c.BasC, c); c.Xb = Xb
    Gp = Xb[:, :c.hp] @ Xb[:, :c.hp].T; Gm = Xb[:, c.hp:] @ Xb[:, c.hp:].T
    c.Gp, c.Gm = Gp, Gm
    B = Gp + Gm; E = math.sqrt(p) * (Gp - Gm)
    c.B = np.rint(B).astype(np.int64); c.E = np.rint(E).astype(np.int64)
    c.resB = np.abs(B - c.B).max(); c.resE = np.abs(E - c.E).max()
    log(f'Gram B integer residual {c.resB:.2e}; E integer residual {c.resE:.2e}; max|B|={np.abs(c.B).max()}, diag min/max {c.B.diagonal().min()}/{c.B.diagonal().max()}')
    Bi = np.linalg.inv(c.B.astype(float))
    c.Binv = Bi
    Neta = c.E @ Bi
    c.Neta = Neta
    log(f'check E B^-1 E = p B: residual {np.abs(c.E @ Bi @ c.E - p * c.B).max():.2e}')
    half = (np.eye(c.n) + Neta) / 2
    log(f'O_F acts ((1+sqrt p)/2 integral matrix): max frac = {np.abs(half-np.rint(half)).max():.2e}; sqrt p integral: {np.abs(Neta-np.rint(Neta)).max():.2e}')
    c.detB = bareiss_det(c.B)
    c.D = math.exp(math.log(c.detB) / (p - 1))
    c.theta = (p - 1) ** 2 / (4 * p)
    log(f'det B = {c.detB}')
    log(f'D = det^(1/(p-1)) = {c.D:.10f}; theta_p = {c.theta:.6f}')
    # factor det a bit
    dd = c.detB; fac = {}
    for q in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]:
        while dd % q == 0:
            dd //= q; fac[q] = fac.get(q, 0) + 1
    log(f'det B factorisation: {fac} cofactor {dd}')
    c.detfac = fac

def coords(c, VC):
    return (remb(VC, c) @ c.Xb.T) @ c.Binv

def gamma_int(c, x):
    """x: (k,n) integer rows -> (Gamma, mplus, mminus, b, e) exact-ish (b,e exact ints)."""
    x = np.atleast_2d(x).astype(np.int64)
    b = np.einsum('ij,jk,ik->i', x, c.B, x); e = np.einsum('ij,jk,ik->i', x, c.E, x)
    sp = math.sqrt(c.p)
    mp = (b + e / sp) / 2; mmn = (b - e / sp) / 2
    return mp * mmn / c.D, mp, mmn, b, e

def fixed_lattice(c, g, order):
    n = c.n
    Cg = coords(c, act(g, c.BasC))
    A = np.rint(Cg).astype(np.int64)
    assert np.abs(Cg - A).max() < 1e-6
    Mx = A - np.eye(n, dtype=np.int64)
    w = 4096.0
    Bk = np.ascontiguousarray(np.hstack([np.eye(n), w * Mx.astype(float)]))
    U = np.ascontiguousarray(np.eye(n, dtype=np.int64)); mu = np.zeros((n, n)); r = np.zeros(n)
    rc = lat.lll_c(n, 2 * n, Bk.ctypes.data, U.ctypes.data, ctypes.c_double(0.99), mu.ctypes.data, r.ctypes.data)
    img = U @ Mx
    isk = (img == 0).all(1)
    F = U[isk]
    rk = int(isk.sum())
    sat = (np.linalg.matrix_rank(img[~isk].astype(float)) == n - rk) if rk < n else True
    assert (F @ A == F).all()
    GF = F @ c.B @ F.T
    U2 = c_lll(GF.astype(float))
    F = U2 @ F
    return F, bool(sat)

def sub_det(c, F):
    return bareiss_det(F.astype(object).dot(c.B.astype(object)).dot(F.astype(object).T))

def grid_enum(c, F, step=0.1, label='', maxnodes=4_000_000_000, tl=1e9):
    """exhaustive: all w in sublattice (rows F) with Gamma(w) < theta must appear. returns min Gamma etc."""
    p = c.p
    a, b, eps = fund_unit(p); le = math.log(eps)
    ts = np.arange(-le, le + step, step)
    GpF = F @ c.Gp @ F.T; GmF = F @ c.Gm @ F.T
    bound0 = 2 * math.sqrt(c.theta * c.D) * math.cosh(step / 2) * (1 + 1e-9)
    best = (1e99, None); tot = 0; allflags = 0; nodes = 0; nviol = 0; viol = {}
    t1 = time.time(); done = 0
    for t in ts:
        if time.time() - t1 > tl:
            break
        A = math.exp(t) * GpF + math.exp(-t) * GmF
        Y, fl, nd = enum_short(A, bound0, maxnodes=maxnodes)
        allflags |= fl; nodes += nd; done += 1
        if len(Y):
            x = Y @ F
            gam, mp, mmn, bb, ee = gamma_int(c, x)
            tot += len(Y)
            i = int(np.argmin(gam))
            if gam[i] < best[0]:
                best = (float(gam[i]), x[i].copy(), float(mp[i]), float(mmn[i]), int(bb[i]), int(ee[i]))
            for j in np.where(gam < c.theta * (1 - 1e-9))[0]:
                viol[tuple(x[j])] = float(gam[j])
    log(f'[{label}] rank {len(F)}: grid {done}/{len(ts)} t-values in [-{le:.4f},{le:.4f}], bound Q_t<= {bound0:.4f}, vectors enumerated {tot}, nodes {nodes}, flags {allflags}, '
        f'min Gamma = {best[0]:.8f} (theta={c.theta:.8f}, ratio {best[0]/c.theta:.6f}), violators (distinct found) {len(viol)}')
    return dict(label=label, rank=int(len(F)), min_gamma=best[0], best=best, complete=(done == len(ts) and allflags == 0),
                n_enum=tot, n_viol=len(viol), viol=viol)

def cert(c, best, extra=None):
    gam, x, mp, mmn, bb, ee = best
    # exact: N = p b^2 - e^2 ; Gamma = N / (4 p D)
    d = dict(p=c.p, x=[int(z) for z in x], B_norm=bb, E_norm=ee, m_plus=mp, m_minus=mmn, D=c.D, Gamma=gam, theta=c.theta,
             N_exact=c.p * bb * bb - ee * ee, formula='Gamma = (p*b^2 - e^2)/(4 p D), b = x B x^T, e = x E x^T, D = det(B)^(1/(p-1))')
    if extra: d.update(extra)
    return d
