"""Equivariant sticky random walk (discrete analogue of Klartag's process) + greedy LP descent
on the G-invariant Ryshkov polyhedron of the fixed module L = O^k (Z-basis coordinates).
Forms: A in Q = {symmetric A : g^T A g = A for g in G}; min_{x in Z^n\\0} A[x] = 1; V = kappa_n/sqrt(det A).
"""
import numpy as np, ctypes, os, json, sys, time, math
from scipy.optimize import linprog

_here = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(_here, 'fp_enum.so'))
_lib.fp_enum.restype = ctypes.c_long
_lib.fp_enum.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_double, ctypes.c_long,
                         ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long,
                         ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_long)]
STATS = {'enum': 0, 'nodes': 0}


def kappa(n):
    return math.pi ** (n / 2) / math.gamma(n / 2 + 1)


def lll_gram(G0, delta=0.99):
    n = G0.shape[0]
    G = np.array(G0, dtype=float)
    U = np.eye(n, dtype=np.int64)
    L = np.linalg.cholesky(G)
    k = 1
    it = 0
    while k < n and it < 100000:
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
            G = U.T @ G0 @ U
            G = (G + G.T) / 2
            L = np.linalg.cholesky(G)
            k = max(k - 1, 1)
    G = U.T @ G0 @ U
    return (G + G.T) / 2, U


def quad(X, A):
    Xf = X.astype(float)
    return np.einsum('ij,ij->i', Xf @ A, Xf)


def sym_basis(n):
    N = n * (n + 1) // 2
    Q = np.zeros((N, n, n))
    k = 0
    for i in range(n):
        for j in range(i, n):
            if i == j:
                Q[k, i, i] = 1.0
            else:
                Q[k, i, j] = Q[k, j, i] = 1 / math.sqrt(2)
            k += 1
    return Q


def qL(u):
    a, b, c, d = u
    return np.array([[a, -b, -c, -d], [b, a, -d, c], [c, d, a, -b], [d, -c, b, a]], dtype=float)


def qmul(p, q):
    return qL(p) @ np.asarray(q, dtype=float)


def cyclo_companion(m):
    import sympy
    x = sympy.symbols('x')
    co = sympy.Poly(sympy.cyclotomic_poly(m, x), x).all_coeffs()[::-1]  # a_0..a_d
    d = len(co) - 1
    M = np.zeros((d, d), dtype=np.int64)
    for j in range(d - 1):
        M[j + 1, j] = 1
    for i in range(d):
        M[i, d - 1] = -int(co[i])
    return M


def icosian_setup():
    s5 = math.sqrt(5)
    tau, sig = (1 + s5) / 2, (1 - s5) / 2
    s = (np.array([.5, .5, .5, .5]), np.array([.5, .5, .5, .5]))
    t = (np.array([tau, sig, 1, 0]) / 2, np.array([sig, tau, 1, 0]) / 2)
    one = (np.array([1., 0, 0, 0]), np.array([1., 0, 0, 0]))
    key = lambda u: tuple(np.round(np.concatenate(u), 8) + 0.0)
    grp = {key(one): one}
    front = [one]
    while front:
        new = []
        for u in front:
            for g in (s, t):
                v = (qmul(g[0], u[0]), qmul(g[1], u[1]))
                kv = key(v)
                if kv not in grp:
                    grp[kv] = v
                    new.append(v)
        front = new
    units = np.array([np.concatenate(u) for u in grp.values()])  # 120 x 8
    assert len(units) == 120, len(units)
    # greedy Z-basis among the units
    idx = []
    for i in range(len(units)):
        if np.linalg.matrix_rank(units[idx + [i]]) == len(idx) + 1:
            idx.append(i)
        if len(idx) == 8:
            break
    Bm = units[idx].T.copy()  # columns = basis
    for it in range(1000):
        C = np.linalg.solve(Bm, units.T)  # 8 x 120 coordinates
        frac = np.abs(C - np.rint(C))
        if frac.max() < 1e-7:
            break
        i, j = np.unravel_index(np.argmax((frac > 1e-7) & (np.abs(C) < 1 - 1e-7)), C.shape)
        if not (frac[i, j] > 1e-7 and abs(C[i, j]) < 1 - 1e-7):
            # reduce: replace using fractional part
            i, j = np.unravel_index(np.argmax(frac), C.shape)
            v = units[j] - Bm @ np.floor(C[:, j])
            Bm[:, i] = v
        else:
            Bm[:, i] = units[j]
    C = np.linalg.solve(Bm, units.T)
    assert np.abs(C - np.rint(C)).max() < 1e-7
    gens = []
    for g in (s, t):
        Lg = np.zeros((8, 8))
        Lg[:4, :4] = qL(g[0])
        Lg[4:, 4:] = qL(g[1])
        M = np.linalg.solve(Bm, Lg @ Bm)
        assert np.abs(M - np.rint(M)).max() < 1e-7
        gens.append(np.rint(M).astype(np.int64))
    return gens, Bm


def family(name):
    """returns (list of generator matrices on O in Z-basis coords, w=|G|, info)"""
    if name == 'eis':
        return [np.array([[0, -1], [1, -1]])], 6, dict(beta=[(1, 0), (-0.5, math.sqrt(3) / 2)])
    if name == 'gauss':
        return [np.array([[0, -1], [1, 0]])], 4, dict(beta=[(1, 0), (0, 1)])
    if name == 'hur':
        Bm = np.array([[1, 0, 0, .5], [0, 1, 0, .5], [0, 0, 1, .5], [0, 0, 0, .5]])
        gens = []
        for u in ([0, 1, 0, 0], [0, 0, 1, 0], [.5, .5, .5, .5]):
            M = np.linalg.solve(Bm, qL(u) @ Bm)
            assert np.abs(M - np.rint(M)).max() < 1e-9
            gens.append(np.rint(M).astype(np.int64))
        return gens, 24, dict(beta=[tuple(Bm[:, i]) for i in range(4)])
    if name.startswith('zeta'):
        m = int(name[4:])
        w = m if m % 2 == 0 else 2 * m
        return [cyclo_companion(m)], w, {}
    if name == 'ico':
        gens, Bm = icosian_setup()
        return gens, 120, dict(Bm=Bm.tolist())
    raise ValueError(name)


class Eq:
    def __init__(self, name, k):
        gens0, w, info = family(name)
        self.name, self.k, self.w, self.info = name, k, w, info
        n0 = gens0[0].shape[0]
        self.n0 = n0
        n = self.n = n0 * k
        self.gens = [np.kron(np.eye(k, dtype=np.int64), g) for g in gens0]
        # group closure
        I = np.eye(n, dtype=np.int64)
        grp = {I.tobytes(): I}
        front = [I]
        while front:
            new = []
            for u in front:
                for g in self.gens:
                    v = g @ u
                    if v.tobytes() not in grp:
                        grp[v.tobytes()] = v
                        new.append(v)
            front = new
        self.group = list(grp.values())
        E = sym_basis(n)
        N = E.shape[0]
        Ef = E.reshape(N, n * n)
        blocks = []
        for g in self.gens:
            gf = g.astype(float)
            T = np.einsum('ji,bjk,kl->bil', gf, E, gf).reshape(N, n * n) - Ef
            blocks.append(T)
        Mst = np.hstack(blocks)
        U, s, Vt = np.linalg.svd(Mst, full_matrices=True)
        r = int((s > 1e-9 * max(s[0], 1)).sum())
        C = U[:, r:].T
        self.QB = (C @ Ef).reshape(-1, n, n)
        self.QBf = self.QB.reshape(-1, n * n)
        self.D = self.QB.shape[0]
        self.U = np.eye(n, dtype=np.int64)
        self.kap = kappa(n)

    # ---------- enumeration
    def short(self, A, bound, maxc=2000, maxnodes=200_000_000, cold=False):
        n = self.n
        if cold:
            G, U = lll_gram(A)
        else:
            G0 = self.U.T @ A @ self.U
            G0 = (G0 + G0.T) / 2
            G, U1 = lll_gram(G0)
            U = self.U @ U1
            self.U = U
        L = np.ascontiguousarray(np.linalg.cholesky(G))
        out = np.zeros((maxc, n), dtype=np.int32)
        vals = np.zeros(maxc)
        flag = ctypes.c_int(0)
        nodes = ctypes.c_long(0)
        cnt = _lib.fp_enum(n, L.ctypes.data, float(bound), maxc, out.ctypes.data, vals.ctypes.data,
                           maxnodes, ctypes.byref(flag), ctypes.byref(nodes))
        STATS['enum'] += 1
        STATS['nodes'] += nodes.value
        if flag.value == 1:
            raise RuntimeError('enumeration node limit')
        Y = out[:cnt].astype(np.int64)
        X = Y @ U.T
        v = quad(X, A) if cnt else np.zeros(0)
        return X, v, flag.value

    def mat(self, c):
        return (c @ self.QBf).reshape(self.n, self.n)

    def phi(self, S):
        Sf = S.astype(float)
        m = len(S)
        out = np.empty((m, self.D))
        for a in range(0, m, 20000):
            Xc = Sf[a:a + 20000]
            out[a:a + 20000] = np.einsum('mi,dij,mj->md', Xc, self.QB, Xc, optimize=True)
        return out

    def null(self, Phi):
        nr = np.linalg.norm(Phi, axis=1)
        P = Phi / nr[:, None]
        if len(P) > 4 * self.D:
            _, ix = np.unique(np.round(P, 6), axis=0, return_index=True)
            P = P[np.sort(ix)]
        u, s, vt = np.linalg.svd(P, full_matrices=True)
        r = int((s > 1e-7 * s[0]).sum())
        return r, vt[r:]

    def polish(self, c, Phi):
        d = np.linalg.lstsq(Phi, 1.0 - Phi @ c, rcond=None)[0]
        return c + d

    def random_start(self, rng):
        n = self.n
        R = rng.standard_normal((n, n))
        M = R @ R.T
        A0 = sum(g.T @ M @ g for g in self.group) / len(self.group)
        c = self.QBf @ A0.ravel()
        A = self.mat(c)
        self.U = np.eye(n, dtype=np.int64)
        G, U = lll_gram(A)
        self.U = U
        b = G.diagonal().min()
        X, v, fl = self.short(A, b * (1 + 1e-9), maxc=20000)
        mn = v.min()
        c = c / mn
        A = A / mn
        S = X[v <= mn * (1 + 1e-9)]
        c = self.polish(c, self.phi(S))
        return c, S

    def line_step(self, c, b, eps=1e-9):
        A = self.mat(c)
        B = self.mat(b)
        Lc = np.linalg.cholesky(A)
        Li = np.linalg.inv(Lc)
        w = np.linalg.eigvalsh(Li @ B @ Li.T)
        if w[0] >= -1e-12:
            raise ValueError('B is PSD')
        tmax = -1.0 / w[0]
        X = None
        for k in range(1, 40):
            tt = tmax * (1 - 0.5 ** k)
            X, v, fl = self.short(A + tt * B, 1 - eps, maxc=2000)
            if len(X):
                break
        if X is None or len(X) == 0:
            raise RuntimeError('no violator found')
        t = tt
        for it in range(100000):
            ax = quad(X, A)
            bx = quad(X, B)
            th = (ax - 1.0) / (-bx)
            t = max(float(th.min()), 0.0)
            X, v, fl = self.short(A + t * B, 1 - eps, maxc=2000)
            if len(X) == 0:
                break
        c2 = c + t * b
        S2, v, fl = self.short(self.mat(c2), 1 + 1e-8, maxc=120000)
        if fl == 2:
            raise RuntimeError('too many minimal vectors')
        c2 = self.polish(c2, self.phi(S2))
        return c2, S2, t

    def to_perfect(self, c, S, rng, mode='random', traj=None):
        steps = 0
        while True:
            r, Fb = self.null(self.phi(S))
            if r >= self.D:
                break
            g = rng.standard_normal(Fb.shape[0])
            b = g @ Fb
            b /= np.linalg.norm(b)
            A = self.mat(c)
            B = self.mat(b)
            if mode == 'descent':
                if np.trace(np.linalg.solve(A, B)) > 0:
                    b = -b
            else:
                Li = np.linalg.inv(np.linalg.cholesky(A))
                if np.linalg.eigvalsh(Li @ B @ Li.T)[0] >= -1e-12:
                    b = -b
            m0 = len(S)
            c, S, t = self.line_step(c, b)
            steps += 1
            if len(S) <= m0 or steps > 20 * self.D + 50:
                raise RuntimeError('walk stuck: no new contact (t=%g)' % t)
            if traj is not None:
                traj.append((r + 1, len(S), self.V(c)))
        return c, S, steps

    def V(self, c):
        sgn, ld = np.linalg.slogdet(self.mat(c))
        return self.kap * math.exp(-0.5 * ld)

    def orbit_rows(self, S):
        Phi = self.phi(S)
        _, ix = np.unique(np.round(Phi, 6), axis=0, return_index=True)
        return Phi[np.sort(ix)]

    def descent_lp(self, c, S):
        A = self.mat(c)
        Pu = self.orbit_rows(S)
        g = self.QBf @ np.linalg.inv(A).ravel()
        res = linprog(g, A_ub=-Pu, b_ub=np.zeros(len(Pu)), A_eq=Pu.sum(0)[None, :], b_eq=[1.0],
                      bounds=[(None, None)] * self.D, method='highs-ds')
        if res.status != 0:
            raise RuntimeError('LP failed: ' + res.message)
        b = res.x
        act = Pu @ b < 1e-7
        Pa = Pu[act]
        d = np.linalg.lstsq(Pa, -Pa @ b, rcond=None)[0]
        b = b + d
        r = np.linalg.matrix_rank(Pa / np.linalg.norm(Pa, axis=1)[:, None], tol=1e-7)
        return float(g @ b), b, int(r), len(Pu)

    def greedy(self, c, S, rng, maxit=200):
        nst = 0
        nonedge = 0
        val = float('nan')
        while nst < maxit:
            val, b, r, norb = self.descent_lp(c, S)
            if val > -1e-9:
                break
            V0 = self.V(c)
            try:
                c2, S2, t = self.line_step(c, b)
            except ValueError:
                break
            if r != self.D - 1:
                nonedge += 1
            if self.null(self.phi(S2))[0] < self.D:
                c2, S2, _ = self.to_perfect(c2, S2, rng, 'descent')
            if not (self.V(c2) > V0 * (1 + 1e-10)):
                break
            c, S = c2, S2
            nst += 1
        return c, S, nst, nonedge, val

    def verify(self, c):
        """independent check: cold LLL, min of the form and number of minimal pairs"""
        A = self.mat(c)
        X, v, fl = self.short(A, 1 - 1e-7, maxc=100, cold=True)
        nviol = len(X)
        X, v, fl = self.short(A, 1 + 1e-7, maxc=120000, cold=True)
        return nviol, len(X), float(v.min()) if len(v) else float('nan'), float(v.max()) if len(v) else float('nan')


def denominator(A, qmax=5000, tol=1e-6):
    for q in range(1, qmax + 1):
        Z = q * A
        if np.abs(Z - np.rint(Z)).max() < tol * q:
            return q
    return 0


def one_run(eq, rng, mode='random'):
    t0 = time.time()
    c, S = eq.random_start(rng)
    V0 = eq.V(c)
    traj = []
    c, S, steps = eq.to_perfect(c, S, rng, mode, traj)
    hw = eq.w // 2
    nv, npairs, vmin, vmax = eq.verify(c)
    rec = dict(V_start=V0, V_perfect=eq.V(c), pairs_perfect=len(S), orbits_perfect=len(S) / hw, steps=steps,
               verify_perfect=[nv, npairs, vmin, vmax],
               traj=[[a, b, round(v, 5)] for a, b, v in traj])
    A = eq.mat(c)
    rec['denom_perfect'] = denominator(A)
    cp = c.copy()
    try:
        c, S, nst, nonedge, val = eq.greedy(c, S, rng)
        nv, npairs, vmin, vmax = eq.verify(c)
        A = eq.mat(c)
        rec.update(V_ext=eq.V(c), pairs_ext=len(S), orbits_ext=len(S) / hw, descent_steps=nst, nonedge=nonedge,
                   lp_val=val, verify_ext=[nv, npairs, vmin, vmax], denom_ext=denominator(A))
    except Exception as e:
        rec['descent_error'] = repr(e)
        c = cp
        A = eq.mat(c)
    rec['A_final'] = np.round(A, 9).tolist()
    rec['time'] = time.time() - t0
    return rec


def run_family(name, k, nruns, tlimit, mode='random', seed=1, outdir=None):
    eq = Eq(name, k)
    rng = np.random.default_rng(abs(hash((name, k, seed, mode))) % (2 ** 32) if False else
                                (sum(ord(ch) for ch in name) * 1000003 + k * 1009 + seed * 17 + (mode == 'descent')))
    recs = []
    fails = 0
    t0 = time.time()
    for i in range(nruns):
        if time.time() - t0 > tlimit:
            break
        try:
            recs.append(one_run(eq, rng, mode))
        except Exception as e:
            fails += 1
            print('FAIL', name, k, i, repr(e), flush=True)
    res = dict(name=name, k=k, n=eq.n, w=eq.w, D=eq.D, mode=mode, seed=seed, fails=fails, elapsed=time.time() - t0,
               stats=dict(STATS), recs=recs)
    outdir = outdir or os.path.join(_here, 'data')
    json.dump(res, open(os.path.join(outdir, f'{name}_k{k}_{mode}_s{seed}.json'), 'w'))
    if recs:
        Vp = np.array([r['V_perfect'] for r in recs])
        Ve = np.array([r.get('V_ext', float('nan')) for r in recs])
        print(f'{name} k={k} n={eq.n} D={eq.D} w={eq.w} mode={mode} runs={len(recs)} fails={fails} '
              f't={time.time()-t0:.1f}s Vperf med={np.median(Vp):.4f} max={Vp.max():.4f} | '
              f'Vext med={np.nanmedian(Ve):.4f} max={np.nanmax(Ve):.4f} '
              f'orb_perf mean={np.mean([r["orbits_perfect"] for r in recs]):.1f}', flush=True)
    else:
        print(f'{name} k={k} n={eq.n} D={eq.D}: no runs, fails={fails}', flush=True)
    return res


if __name__ == '__main__':
    # usage: eq.py mode seed name:k:nruns:tlimit ...
    # usage: eq.py seed mode:name:k:nruns:tlimit ...
    seed = int(sys.argv[1])
    for spec in sys.argv[2:]:
        mode, name, k, nr, tl = spec.split(':')
        try:
            run_family(name, int(k), int(nr), float(tl), mode, seed)
        except Exception as e:
            print('FAMILY FAIL', spec, repr(e), flush=True)
