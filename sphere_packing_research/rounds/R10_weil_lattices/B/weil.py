#!/usr/bin/env python3
# Builder B: exact-integer construction of the Weil lattice Lambda_0 = Z[SL_2(p)] v0 (p = 5 mod 8) and Gamma tests.
import numpy as np, ctypes, math, sys, json, time, os
HERE = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(os.path.join(HERE, 'lat.so'))
I64 = np.ctypeslib.ndpointer(np.int64, flags='C_CONTIGUOUS')
lib.lll_gram.argtypes = [ctypes.c_int, I64, I64, ctypes.c_double]; lib.lll_gram.restype = ctypes.c_long
lib.fp_enum.argtypes = [ctypes.c_int, np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS'), ctypes.c_double, ctypes.c_long,
                        np.ctypeslib.ndpointer(np.int32, flags='C_CONTIGUOUS'), np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS'),
                        ctypes.c_long, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_long)]
lib.fp_enum.restype = ctypes.c_long
T0 = time.time()
def log(*a):
    print('[%7.1fs]' % (time.time() - T0), *a, flush=True)

def lll_int(G, delta=0.99):
    n = G.shape[0]
    Gc = np.ascontiguousarray(np.array(G, dtype=np.int64)); U = np.ascontiguousarray(np.eye(n, dtype=np.int64))
    lib.lll_gram(n, Gc, U, delta)
    return U, Gc

def fp(Q, bound, maxc=400000, maxnodes=4*10**9):
    n = Q.shape[0]
    Lc = np.ascontiguousarray(np.linalg.cholesky(Q))
    out = np.zeros((maxc, n), dtype=np.int32); vals = np.zeros(maxc)
    flag = ctypes.c_int(0); nodes = ctypes.c_long(0)
    cnt = lib.fp_enum(n, Lc, bound, maxc, out, vals, maxnodes, ctypes.byref(flag), ctypes.byref(nodes))
    return out[:cnt].astype(np.int64), vals[:cnt], flag.value, nodes.value

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def prim_root(p):
    fs = [q for q in range(2, p) if (p - 1) % q == 0 and all(q % d for d in range(2, q))]
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fs): return g

def fund_unit(p):
    y = 1
    while True:
        for s in (-4, 4):
            x2 = p * y * y + s; x = math.isqrt(x2)
            if x * x == x2: return (x + y * math.sqrt(p)) / 2, (x, y, s)
        y += 1

def det_exact(M):
    M = [[int(x) for x in row] for row in M]; n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    M[k], M[i] = M[i], M[k]; sign = -sign; break
            else: return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def hnf_mod(R, Dn):
    """upper triangular basis (n x n) of rowspace_Z(R) + Dn Z^n"""
    R = np.array(R, dtype=np.int64) % Dn
    n = R.shape[1]
    R = R[np.any(R != 0, axis=1)]
    out = np.zeros((n, n), dtype=np.int64)
    for j in range(n):
        piv = np.zeros((1, n), dtype=np.int64); piv[0, j] = Dn
        nz = R[:, j] != 0
        Wk = np.vstack([R[nz], piv]); rest = R[~nz]
        while True:
            col = Wk[:, j]; nzi = np.nonzero(col)[0]
            if len(nzi) == 1: break
            i0 = nzi[np.argmin(col[nzi])]; g = col[i0]
            q = col // g; q[i0] = 0
            Wk = Wk - q[:, None] * Wk[i0][None, :]
            Wk[:, j + 1:] %= Dn
        i0 = np.nonzero(Wk[:, j])[0][0]
        out[j] = Wk[i0]
        others = np.delete(Wk, i0, axis=0)
        others = others[np.any(others != 0, axis=1)]
        R = np.vstack([rest, others])
    return out

def exdiv(M, d):
    assert not np.any(M % d), 'inexact division'
    return M // d

class Weil:
    def __init__(self, p):
        assert p % 8 == 5
        self.p = p; h = (p - 1) // 2; m = (p - 1) // 4; n = 2 * (p - 1); Dn = m * p
        self.h, self.m, self.n, self.Dn = h, m, n, Dn
        g = prim_root(p); self.g = g
        mu = sorted({pow(g, 4 * i, p) for i in range(m)})
        zt = np.exp(2j * np.pi * np.arange(p) / p)
        a = np.arange(1, h + 1); x = np.arange(1, h + 1)
        chi = np.array([legendre(int(aa), p) for aa in a]); sp = math.sqrt(p)
        ex = (2 * a[:, None, None] * x[None, :, None] * x[None, None, :]) % p
        S0 = (zt[ex] - zt[(-ex) % p]) / (chi[:, None, None] * sp)
        Td = zt[(a[:, None] * (x[None, :] ** 2)) % p]
        M = S0[0] @ np.diag(Td[0]); M3 = M @ M @ M
        if np.allclose(M3, -np.eye(h)): kappa = 1
        elif np.allclose(M3, np.eye(h)): kappa = -1
        else: raise RuntimeError('(ST)^3 not +-I')
        self.kappa = kappa
        S_all = kappa * S0
        for i in range(h):
            M = S_all[i] @ np.diag(Td[i])
            assert np.allclose(M @ M @ M, -np.eye(h)) and np.allclose(S_all[i] @ S_all[i], -np.eye(h))
        W = np.zeros((2, h))
        for ci, c in enumerate([1, g]):
            for u in mu:
                y = (c * u) % p
                if y <= h: W[ci, y - 1] += 1
                else: W[ci, p - y - 1] -= 1
        Bas = np.zeros((n, h, h), complex)
        for ci in range(2):
            for jj in range(1, p):
                Bas[ci * (p - 1) + jj - 1] = W[ci][None, :] * zt[(a[:, None] * jj * (x[None, :] ** 2)) % p]
        def ip(U, V): return 2 * np.real(U.reshape(len(U), -1).conj() @ V.reshape(len(V), -1).T)
        def rnd(Mf, name):
            Mi = np.rint(Mf); res = np.abs(Mf - Mi).max()
            assert res < 1e-6, (name, res)
            return Mi.astype(np.int64), res
        GB, r1 = rnd(ip(Bas, Bas), 'GB')
        J = np.ones((p - 1, p - 1), dtype=np.int64); I = np.eye(p - 1, dtype=np.int64); Z = np.zeros((p - 1, p - 1), dtype=np.int64)
        GBx = np.block([[m * (p * I - J), Z], [Z, m * (p * I - J)]])
        assert np.array_equal(GB, GBx)
        IJ = np.block([[I + J, Z], [Z, I + J]])
        etaBas = Bas * (chi[:, None] * sp)[None, :, :]
        EB, r2 = rnd(ip(etaBas, Bas), 'EB')
        assert np.array_equal(EB, EB.T)
        SBas = np.einsum('axy,iay->iax', S_all, Bas)
        YS, r3 = rnd(ip(SBas, Bas), 'YS')
        NS = YS @ IJ
        Meta = exdiv(EB @ IJ, Dn)
        MT = np.zeros((n, n), dtype=np.int64)
        for ci in range(2):
            o = ci * (p - 1)
            for jj in range(1, p - 1): MT[o + jj - 1, o + jj] = 1
            MT[o + p - 2, o:o + p - 1] = -1
        self.round_res = (r1, r2, r3)
        self.Bas, self.S_all, self.Td, self.chi, self.Wc = Bas, S_all, Td, chi, W
        # exact relation checks (python ints)
        O = lambda A: A.astype(object)
        In = np.eye(n, dtype=object)
        NSo, MTo, GBo, Meo = O(NS), O(MT), O(GB), O(Meta)
        assert np.array_equal(NSo @ NSo, -Dn * Dn * In), 'S^2'
        TS = MTo @ NSo
        assert np.array_equal(TS @ TS @ TS, -Dn ** 3 * In), '(ST)^3'
        Tp = In
        for _ in range(p): Tp = Tp @ MTo
        assert np.array_equal(Tp, In), 'T^p'
        assert np.array_equal(NSo @ GBo @ NSo.T, Dn * Dn * GBo), 'S isometry'
        assert np.array_equal(MTo @ GBo @ MTo.T, GBo), 'T isometry'
        assert np.array_equal(Meo @ Meo, p * In), 'eta^2=p'
        assert np.array_equal(Meo @ MTo, MTo @ Meo) and np.array_equal(Meo @ NSo, NSo @ Meo), 'eta central'
        assert np.array_equal(Meo @ GBo, O(EB)), 'E = eta G'
        self.GB, self.EB, self.NS, self.MT, self.Meta = GB, EB, NS, MT, Meta
        self.eps, self.eps_xy = fund_unit(p)
        self.theta = (p - 1) ** 2 / (4 * p)
        log('p=%d setup ok: kappa=%d, rounding residuals %s, exact relations S^2=-1,(ST)^3=-1,T^p=1,isometry,eta^2=p,eta central verified' % (p, kappa, self.round_res))

    def build(self):
        p, n, Dn, m = self.p, self.n, self.Dn, self.m
        L = hnf_mod(self.NS, Dn); it = 0
        while True:
            idx = [int(v) for v in np.diag(L)]
            R = np.vstack([L, L @ self.MT, exdiv(L @ self.NS, Dn)])
            L2 = hnf_mod(R, Dn); it += 1
            if [int(v) for v in np.diag(L2)] == idx: break
            L = L2
        self.A = L
        dprod = 1
        for v in np.diag(L): dprod *= int(v)
        detGB = (m ** (p - 1) * p ** (p - 2)) ** 2
        num = detGB * dprod * dprod; den = Dn ** (2 * n)
        assert num % den == 0
        self.detB = num // den
        self.Dexact = None
        if self.detB == (m * m * p) ** (p - 1): self.Dexact = m * m * p
        self.D = float(self.Dexact) if self.Dexact else math.exp(math.log(self.detB) / (p - 1))
        Bg = exdiv(L @ self.GB @ L.T, Dn * Dn); Eg = exdiv(L @ self.EB @ L.T, Dn * Dn)
        log('closure iterations %d; det(B on Lambda_0) = %s ; D = %s (exact m^2 p: %s)' % (it, self.detB if self.detB < 10**40 else '~1e%d' % (len(str(self.detB)) - 1), self.D, self.Dexact is not None))
        U, Br = lll_int(Bg, 0.99)
        Ar = U @ L
        assert np.array_equal(exdiv(Ar @ self.GB @ Ar.T, Dn * Dn), Br)
        Er = exdiv(Ar @ self.EB @ Ar.T, Dn * Dn)
        self.Ar, self.B, self.E = Ar, Br, Er
        # sanity E B^-1 E = p B (float)
        Bf = Br.astype(float); Ef = Er.astype(float)
        res = np.abs(Ef @ np.linalg.solve(Bf, Ef) - p * Bf).max()
        log('LLL done; diag(B) min %d max %d ; |E B^-1 E - pB| = %.2e' % (Br.diagonal().min(), Br.diagonal().max(), res))
        return self

    def isometry(self, X):
        """X = images (scaled coords) of the basis rows Ar; returns integer P with P Ar = X (verified exactly)"""
        Pf = np.linalg.solve(self.Ar.T.astype(float), X.T.astype(float)).T
        P = np.rint(Pf).astype(np.int64)
        assert np.abs(Pf - P).max() < 1e-4 and np.array_equal(P @ self.Ar, X), 'not an integral endomorphism of the lattice'
        return P

    def torus_perm(self, q):
        """isometry P of generator of the order-q subgroup of the split torus H (permutation of basis b_{j,c})"""
        p, m = self.p, self.m
        a0 = pow(self.g, 4 * (m // q), p); s = pow(pow(a0, 2, p), p - 2, p)
        X = np.zeros_like(self.Ar)
        for ci in range(2):
            o = ci * (p - 1)
            for jj in range(1, p):
                X[:, o + (jj * s) % p - 1] = self.Ar[:, o + jj - 1]
        return self.isometry(X)

    def word(self, w):
        """isometry of a word in S,T (applied left to right on row coords)"""
        X = self.Ar.copy()
        for ch in w:
            X = X @ self.MT if ch == 'T' else exdiv(X @ self.NS, self.Dn)
        return self.isometry(X)

    def fixed(self, P, lam=10**5):
        n = self.n; K = P - np.eye(n, dtype=np.int64)
        if not np.any(K): return np.eye(n, dtype=np.int64)
        Gk = np.eye(n, dtype=np.int64) + lam * (K @ K.T)
        U, _ = lll_int(Gk, 0.99)
        ker = np.array([not np.any(U[i] @ K) for i in range(n)])
        k = int(ker.sum())
        assert ker[:k].all() and not ker[k:].any(), 'kernel vectors not leading'
        C = U[:k]
        Bs = C @ self.B @ C.T
        U2, _ = lll_int(Bs, 0.99)
        return U2 @ C

    def fixed_perm(self, q):
        """EXACT fixed sublattice of the order-q subgroup of the split torus (permutation action), via HNF mod Dn in difference coordinates"""
        p, m, n, Dn = self.p, self.m, self.n, self.Dn
        a0 = pow(self.g, 4 * (m // q), p); s = pow(a0, 2, p)
        seen = set(); orbs = []
        for ci in range(2):
            o = ci * (p - 1)
            for jj in range(1, p):
                if (ci, jj) in seen: continue
                orb = []; x = jj
                while (ci, x) not in seen:
                    seen.add((ci, x)); orb.append(o + x - 1); x = (x * s) % p
                orbs.append(orb)
        k = len(orbs); A = self.A
        diffs = []; reps = []
        for orb in orbs:
            reps.append(A[:, orb[0]])
            for j in orb[1:]: diffs.append(A[:, j] - A[:, orb[0]])
        Lp = np.column_stack(diffs + reps)
        Hn = hnf_mod(Lp, Dn)
        assert not np.any(Hn[n - k:, :n - k])
        tpart = Hn[n - k:, n - k:]
        X = np.zeros((k, n), dtype=np.int64)
        for i, orb in enumerate(orbs): X[:, orb] = tpart[:, [i]]
        C = self.isometry(X)
        U2, _ = lll_int(C @ self.B @ C.T, 0.99)
        return U2 @ C

    def grams(self, C):
        return C @ self.B @ C.T, C @ self.E @ C.T

    def gamma_data(self, Wv, Bs, Es):
        b = np.einsum('ij,jk,ik->i', Wv, Bs, Wv); e = np.einsum('ij,jk,ik->i', Wv, Es, Wv)
        N = self.p * b * b - e * e          # Gamma = N / (4 D p)
        return b, e, N

    def tsearch(self, C, step=0.1, do_enum=True, maxc=400000, maxnodes=4 * 10**9, label='', tlim=None, pairs=False, trange=None):
        """C: k x n integer basis (coords in reduced basis) of a sublattice. Exhaustive (if do_enum) search for Gamma<theta."""
        p = self.p; sp = math.sqrt(p); D = self.D; theta = self.theta
        Bs, Es = self.grams(C); k = len(C)
        le = math.log(self.eps); nt = int(math.ceil(2 * le / step)); ts = [-le + step * (i + 0.5) for i in range(nt)]
        bound = 2 * math.sqrt(theta * D) * math.cosh(step / 2) * (1 + 1e-9)
        thrN = D * (p - 1) ** 2
        U = np.eye(k, dtype=np.int64)
        best = None; nviol = 0; nborder = 0; complete = True; totnodes = 0; tstart = time.time(); done_t = 0
        viol = None
        if trange: ts = ts[trange[0]:trange[1]]
        for t in ts:
            if tlim and time.time() - tstart > tlim: complete = False; break
            def Qof(U):
                Bc = (U @ Bs @ U.T).astype(float); Ec = (U @ Es @ U.T).astype(float)
                return math.exp(t) * (Bc + Ec / sp) / 2 + math.exp(-t) * (Bc - Ec / sp) / 2
            Q = Qof(U)
            Qi = np.rint(Q * (2.0 ** 36 / np.abs(Q).max())).astype(np.int64)
            U2, _ = lll_int(Qi, 0.99); U = U2 @ U; Q = Qof(U)
            cands = [U]
            if pairs:
                iu = np.triu_indices(k, 1)
                cands.append(U[iu[0]] + U[iu[1]]); cands.append(U[iu[0]] - U[iu[1]])
            if do_enum:
                X, vals, flag, nodes = fp(Q, bound, maxc, maxnodes); totnodes += nodes
                if flag: complete = False
                if len(X): cands.append(X @ U)
            Wv = np.vstack(cands)
            b, e, N = self.gamma_data(Wv, Bs, Es)
            i = int(np.argmin(N)); gam = N[i] / (4 * D * p)
            if best is None or gam < best[0]: best = (float(gam), int(b[i]), int(e[i]), t, (Wv[i] @ C).tolist())
            nv = int((N < thrN * (1 - 1e-12)).sum()); nviol += nv; nborder += int((np.abs(N - thrN) <= thrN * 1e-12).sum())
            done_t += 1
            if trange or k >= 100: log('   t=%.3f done: best so far Gamma/theta=%.5f, violator hits %d, complete=%s, nodes %.2e' % (t, best[0] / theta, nviol, complete, totnodes))
        res = dict(label=label, rank=k, t_points=nt, t_done=done_t, complete=complete, bound=bound, min_gamma=best[0], theta=theta,
                   ratio=best[0] / theta, b=best[1], e=best[2], t_at_min=best[3], n_violator_hits=nviol, n_border_hits=nborder, nodes=totnodes, w=best[4])
        log('%s rank %d: t-grid %d/%d complete=%s minGamma=%.6f theta=%.6f ratio=%.6f (b=%d,e=%d) violator hits=%d border hits=%d nodes=%.2e' %
            (label, k, done_t, nt, complete, best[0], theta, best[0] / theta, best[1], best[2], nviol, nborder, totnodes))
        return res

def run(p, tasks):
    out = {}
    Wl = Weil(p).build(); m, n = Wl.m, Wl.n
    np.save(os.path.join(HERE, 'B_p%d.npy' % p), Wl.B); np.save(os.path.join(HERE, 'E_p%d.npy' % p), Wl.E)
    np.save(os.path.join(HERE, 'Ared_p%d.npy' % p), Wl.Ar)
    out.update(p=p, n=n, m=m, detB=str(Wl.detB), D=Wl.D, D_exact=Wl.Dexact, theta=Wl.theta, eps=Wl.eps, kappa=Wl.kappa)
    # H-fixed lattice and order
    if m > 1:
        CH = Wl.fixed_perm(m)
    else:
        CH = np.eye(n, dtype=np.int64)
    BH, EH = Wl.grams(CH); dH = det_exact(BH)
    out['rank_H'] = len(CH); out['det_BH'] = str(dH); out['order_maximal'] = (dH == p ** 4 * m ** 16)
    out['calibration_H_borderline'] = (dH ** m == m ** (8 * m) * Wl.detB)
    log('Lambda^H rank %d det_B = %d ; maximal-order value p^4 m^16 = %d -> maximal: %s ; ratio det/(p^4 m^16) = %s ; borderline calibration det^m == m^(8m) det(Lambda): %s' %
        (len(CH), dH, p ** 4 * m ** 16, out['order_maximal'], dH / (p ** 4 * m ** 16), out['calibration_H_borderline']))
    Wl.CH = CH
    return Wl, out

if __name__ == '__main__':
    p = int(sys.argv[1])
    Wl, out = run(p, sys.argv[2:])
    json.dump(out, open(os.path.join(HERE, 'basic_p%d.json' % p), 'w'), indent=1)
