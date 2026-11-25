"""Random perfect forms (sticky random walk, discrete analogue of Klartag's process)
and greedy descent to extreme forms on the Ryshkov polyhedron of Z^n.
Forms A are normalized so that min_{x in Z^n \\ 0} A[x] = 1.  V = kappa_n / sqrt(det A).
"""
import numpy as np, ctypes, os, json, sys, time, math
from scipy.optimize import linprog

_here = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(_here, 'fp_enum.so'))
_lib.fp_enum.restype = ctypes.c_long
_lib.fp_enum.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_double, ctypes.c_long,
                         ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long,
                         ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_long)]
_MAXC = 20000
_out = np.zeros((_MAXC, 16), dtype=np.int32)
_vals = np.zeros(_MAXC)
STATS = {'enum': 0, 'nodes': 0}


def kappa(n):
    return math.pi ** (n / 2) / math.gamma(n / 2 + 1)


def lll_gram(G0, delta=0.99):
    """LLL on a Gram matrix. Returns (G, U) with G = U^T G0 U, U unimodular integer."""
    n = G0.shape[0]
    G = np.array(G0, dtype=float)
    U = np.eye(n, dtype=np.int64)
    L = np.linalg.cholesky(G)
    k = 1
    it = 0
    while k < n and it < 20000:
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


def short_vectors(A, bound, maxc=5000, maxnodes=50_000_000):
    """All nonzero integer x (one per +-pair) with A[x] <= bound.  Returns X (m,n), vals, flag."""
    n = A.shape[0]
    G, U = lll_gram(A)
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
    v = np.einsum('ij,jk,ik->i', X, A, X) if cnt else np.zeros(0)
    return X, v, flag.value


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
    return Q.reshape(N, n * n)


_QB = {}


def contact_matrix(S):
    m, n = S.shape
    return np.einsum('ki,kj->kij', S, S).reshape(m, n * n).astype(float)


def F_basis(S):
    """Orthonormal basis (rows, as n*n vectors) of F = {B sym : B[x]=0 for x in S}; also rank."""
    n = S.shape[1]
    if n not in _QB:
        _QB[n] = sym_basis(n)
    Q = _QB[n]
    Ms = contact_matrix(S) @ Q.T
    u, s, vt = np.linalg.svd(Ms, full_matrices=True)
    r = int((s > 1e-9 * s[0]).sum())
    return vt[r:] @ Q, r


def polish(A, S):
    n = A.shape[0]
    M = contact_matrix(S)
    d = np.linalg.lstsq(M, 1.0 - M @ A.ravel(), rcond=None)[0]
    A2 = A + d.reshape(n, n)
    return (A2 + A2.T) / 2


def rebase(A, S):
    G, U = lll_gram(A)
    Ui = np.rint(np.linalg.inv(U)).astype(np.int64)
    return G, S @ Ui.T


def line_step(A, B, eps=1e-9):
    """Largest t>0 with min(A+tB) >= 1 (requires B not PSD). Returns A_new, S_new, t."""
    Lc = np.linalg.cholesky(A)
    Li = np.linalg.inv(Lc)
    w = np.linalg.eigvalsh(Li @ B @ Li.T)
    if w[0] >= -1e-12:
        raise ValueError('B is PSD')
    tmax = -1.0 / w[0]
    X = None
    for k in range(1, 45):
        tt = tmax * (1 - 0.5 ** k)
        X, v, fl = short_vectors(A + tt * B, 1 - eps, maxc=2000)
        if len(X):
            break
    if X is None or len(X) == 0:
        raise RuntimeError('no violator found')
    t = tt
    for it in range(100000):
        ax = np.einsum('ij,jk,ik->i', X, A, X)
        bx = np.einsum('ij,jk,ik->i', X, B, X)
        th = (ax - 1.0) / (-bx)
        t = max(float(th.min()), 0.0)
        X, v, fl = short_vectors(A + t * B, 1 - eps, maxc=2000)
        if len(X) == 0:
            break
    A2 = A + t * B
    S2, v, fl = short_vectors(A2, 1 + 1e-8, maxc=5000)
    A2 = polish(A2, S2)
    A2, S2 = rebase(A2, S2)
    return A2, S2, t


def to_perfect(A, S, rng, mode='descent', traj=None):
    """Sticky random walk in F until the form is perfect."""
    n = A.shape[0]
    N = n * (n + 1) // 2
    steps = 0
    while True:
        Fb, r = F_basis(S)
        if r == N:
            break
        g = rng.standard_normal(Fb.shape[0])
        B = (g @ Fb).reshape(n, n)
        B = (B + B.T) / 2
        B /= np.linalg.norm(B)
        c = np.trace(np.linalg.solve(A, B))
        if mode == 'descent':
            if c > 0:
                B = -B
        else:  # random sign; flip only if PSD (infinite ray)
            Lc = np.linalg.cholesky(A)
            Li = np.linalg.inv(Lc)
            if np.linalg.eigvalsh(Li @ B @ Li.T)[0] >= -1e-12:
                B = -B
        A, S, t = line_step(A, B)
        steps += 1
        if traj is not None:
            traj.append((int(F_basis(S)[1]), len(S), kappa(n) / math.sqrt(np.linalg.det(A))))
    return A, S, steps


def random_start(n, rng):
    G = rng.standard_normal((n, n))
    A, U = lll_gram(G.T @ G)
    b = A.diagonal().min()
    X, v, fl = short_vectors(A, b * (1 + 1e-12))
    m = v.min()
    A = A / m
    S = X[v <= m * (1 + 1e-9)]
    A = polish(A, S)
    return A, S


def descent_lp(A, S):
    """min Tr(A^-1 B) s.t. B[x]>=0 (x in S), sum_x B[x] = 1.  Value = max-min eutaxy coefficient."""
    n = A.shape[0]
    iu = np.triu_indices(n)
    w = np.where(iu[0] == iu[1], 1.0, 2.0)
    Phi = (S[:, iu[0]] * S[:, iu[1]]) * w
    Ainv = np.linalg.inv(A)
    c = Ainv[iu] * w
    res = linprog(c, A_ub=-Phi, b_ub=np.zeros(len(S)), A_eq=Phi.sum(0)[None, :], b_eq=[1.0],
                  bounds=[(None, None)] * len(c), method='highs-ds')
    if res.status != 0:
        raise RuntimeError('LP failed: ' + res.message)
    B = np.zeros((n, n))
    B[iu] = res.x
    B = B + B.T - np.diag(np.diag(B))
    bx = np.einsum('ij,jk,ik->i', S, B, S)
    act = bx < 1e-7
    M = contact_matrix(S[act])
    d = np.linalg.lstsq(M, -M @ B.ravel(), rcond=None)[0]
    B = B + d.reshape(n, n)
    B = (B + B.T) / 2
    Fb, r = F_basis(S[act])
    val = float(np.trace(Ainv @ B))
    return val, B, r


def eutaxy_lp(A, S):
    """max t s.t. A^-1 = sum lam_x x x^T, lam_x >= t."""
    n = A.shape[0]
    iu = np.triu_indices(n)
    Phi = (S[:, iu[0]] * S[:, iu[1]]).astype(float)
    m = len(S)
    Aeq = np.hstack([Phi.T, np.zeros((len(iu[0]), 1))])
    beq = np.linalg.inv(A)[iu]
    Aub = np.hstack([-np.eye(m), np.ones((m, 1))])
    c = np.zeros(m + 1)
    c[-1] = -1
    res = linprog(c, A_ub=Aub, b_ub=np.zeros(m), A_eq=Aeq, b_eq=beq, bounds=[(None, None)] * (m + 1),
                  method='highs-ds')
    if res.status != 0:
        return float('nan')
    return float(res.x[-1])


def greedy_descent(A, S, rng, maxit=500):
    n = A.shape[0]
    N = n * (n + 1) // 2
    nsteps = 0
    nonedge = 0
    Vpath = [kappa(n) / math.sqrt(np.linalg.det(A))]
    while nsteps < maxit:
        val, B, r = descent_lp(A, S)
        if val > 1e-9:
            break
        d0 = np.linalg.det(A)
        A2, S2, t = line_step(A, B)
        if r != N - 1:
            nonedge += 1
        if F_basis(S2)[1] < N:
            A2, S2, _ = to_perfect(A2, S2, rng, 'descent')
        d1 = np.linalg.det(A2)
        if not (d1 < d0 * (1 - 1e-12)):
            break
        A, S = A2, S2
        nsteps += 1
        Vpath.append(kappa(n) / math.sqrt(d1))
    return A, S, nsteps, nonedge, val, Vpath


def denominator(A, qmax=20000, tol=1e-6):
    for q in range(1, qmax + 1):
        Z = q * A
        if np.abs(Z - np.rint(Z)).max() < tol:
            return q
    return 0


def ip_invariant(A, S):
    """histogram of |x^T A y| over pairs of minimal vectors (finer invariant)."""
    Gm = np.abs(S @ A @ S.T)
    iu = np.triu_indices(len(S), 1)
    vals, cnts = np.unique(np.round(Gm[iu], 5), return_counts=True)
    return [[float(a), int(b)] for a, b in zip(vals, cnts)]


def one_sample(n, rng, mode='descent'):
    t0 = time.time()
    A, S = random_start(n, rng)
    V0 = kappa(n) / math.sqrt(np.linalg.det(A))
    traj = []
    A, S, steps = to_perfect(A, S, rng, mode, traj)
    detp = float(np.linalg.det(A))
    rec = dict(n=n, mode=mode, V_start=V0, V_perfect=kappa(n) / math.sqrt(detp), det_perfect=detp,
               pairs_perfect=len(S), steps_gen=steps, denom_perfect=denominator(A),
               traj=[[a, b, round(c, 6)] for a, b, c in traj])
    lpv, _, _ = descent_lp(A, S)
    rec['eut_margin_perfect'] = lpv
    A, S, nst, nonedge, val, Vpath = greedy_descent(A, S, rng)
    dete = float(np.linalg.det(A))
    rec.update(V_ext=kappa(n) / math.sqrt(dete), det_ext=dete, pairs_ext=len(S), descent_steps=nst,
               nonedge=nonedge, lp_val_ext=val, eut_lp=eutaxy_lp(A, S), denom_ext=denominator(A),
               Vpath=[round(x, 6) for x in Vpath], ip_inv=ip_invariant(A, S) if n <= 9 else None,
               A_ext=np.round(A, 10).tolist(), time=time.time() - t0)
    return rec


if __name__ == '__main__':
    n = int(sys.argv[1])
    nsamp = int(sys.argv[2])
    tlimit = float(sys.argv[3])
    mode = sys.argv[4] if len(sys.argv) > 4 else 'descent'
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    rng = np.random.default_rng(1000 * n + seed)
    recs = []
    t0 = time.time()
    fails = 0
    for i in range(nsamp):
        if time.time() - t0 > tlimit:
            break
        try:
            recs.append(one_sample(n, rng, mode))
        except Exception as e:
            fails += 1
            print('FAIL', n, i, repr(e), flush=True)
    out = os.path.join(_here, 'data', f'runs_n{n}_{mode}_s{seed}.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(dict(n=n, mode=mode, seed=seed, fails=fails, elapsed=time.time() - t0, stats=STATS, recs=recs),
              open(out, 'w'))
    Vp = np.array([r['V_perfect'] for r in recs])
    Ve = np.array([r['V_ext'] for r in recs])
    print(f'n={n} mode={mode} samples={len(recs)} fails={fails} time={time.time()-t0:.1f}s '
          f'Vperf med={np.median(Vp):.4f} max={Vp.max():.4f} | Vext med={np.median(Ve):.4f} max={Ve.max():.4f} '
          f'pairs_ext mean={np.mean([r["pairs_ext"] for r in recs]):.1f}', flush=True)
