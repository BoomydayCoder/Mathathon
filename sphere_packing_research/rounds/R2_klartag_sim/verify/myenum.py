#!/usr/bin/env python3
"""Verifier's own (from-scratch) lattice enumeration: LLL on a Gram matrix (integer transform tracked),
depth-first Fincke-Pohst, and a completely different rigorous brute-force box enumeration
(|z_i| <= sqrt(bound * (G^-1)_ii) by Cauchy-Schwarz)."""
import math, itertools
import numpy as np

def kappa(n): return math.pi ** (n / 2) / math.gamma(n / 2 + 1)

def my_lll(Bm, delta=0.75):
    """Textbook LLL on the columns of Bm (float). Returns integer unimodular U (object ints) with Bm@U reduced."""
    B = np.array(Bm, dtype=float).copy()
    n = B.shape[1]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]   # columns = new basis in old coords
    def gso():
        Bs = np.zeros_like(B); mu = np.zeros((n, n)); nb = np.zeros(n)
        for i in range(n):
            v = B[:, i].copy()
            for j in range(i):
                mu[i, j] = (B[:, i] @ Bs[:, j]) / nb[j]
                v -= mu[i, j] * Bs[:, j]
            Bs[:, i] = v; nb[i] = v @ v
        return mu, nb
    k = 1; it = 0
    mu, nb = gso()
    while k < n and it < 100000:
        it += 1
        for j in range(k - 1, -1, -1):
            q = int(round(mu[k, j]))
            if q != 0:
                B[:, k] -= q * B[:, j]
                for r in range(n): U[r][k] -= q * U[r][j]
                mu, nb = gso()
        if nb[k] >= (delta - mu[k, k - 1] ** 2) * nb[k - 1]:
            k += 1
        else:
            B[:, [k - 1, k]] = B[:, [k, k - 1]]
            for r in range(n): U[r][k - 1], U[r][k] = U[r][k], U[r][k - 1]
            mu, nb = gso()
            k = max(k - 1, 1)
    return U

def int_det(U):
    """exact determinant of an integer matrix (fraction-free Bareiss)."""
    M = [list(map(int, r)) for r in U]; n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None: return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def fp_dfs(G, bound):
    """All nonzero z in Z^n with z^T G z <= bound, one per +-pair. Depth-first Fincke-Pohst (own code)."""
    n = G.shape[0]
    R = np.linalg.cholesky(G).T
    res = []; z = [0] * n
    slack = bound * 1e-9 + 1e-12
    def rec(i, partial):
        if i < 0:
            res.append(tuple(z)); return
        c = 0.0
        for j in range(i + 1, n): c -= R[i, j] * z[j]
        c /= R[i, i]
        rad = math.sqrt(max(bound + slack - partial, 0.0)) / R[i, i]
        lo = math.ceil(c - rad - 1e-9); hi = math.floor(c + rad + 1e-9)
        for zi in range(lo, hi + 1):
            z[i] = zi
            t = R[i, i] * (zi - c)
            p2 = partial + t * t
            if p2 <= bound + slack: rec(i - 1, p2)
        z[i] = 0
    rec(n - 1, 0.0)
    out = []
    for v in res:
        nzv = [a for a in v if a != 0]
        if nzv and nzv[0] > 0: out.append(v)
    Z = np.array(out, dtype=float).reshape(-1, n)
    q = np.einsum('ij,jk,ik->i', Z, G, Z) if len(Z) else np.zeros(0)
    keep = q <= bound + slack
    return Z[keep], q[keep]

def box_enum(G, bound, cap=3e6):
    """Rigorous brute force: z^T G z <= bound implies |z_i| <= sqrt(bound*(G^-1)_ii). Returns (Z,q) up to sign or None if too big."""
    n = G.shape[0]
    Gi = np.linalg.inv(G)
    b = np.floor(np.sqrt(bound * (1 + 1e-9) * np.diag(Gi)) + 1e-9).astype(int)
    size = float(np.prod(2.0 * b + 1))
    if size > cap: return None, b, size
    # split: leading coordinates by itertools, trailing by meshgrid
    order = list(range(n))
    tail = []; prod = 1
    for i in reversed(order):
        if prod * (2 * b[i] + 1) <= 300000: tail.insert(0, i); prod *= (2 * b[i] + 1)
        else: break
    head = [i for i in order if i not in tail]
    grids = np.meshgrid(*[np.arange(-b[i], b[i] + 1) for i in tail], indexing='ij') if tail else []
    T = np.stack([g.ravel() for g in grids], axis=1).astype(float) if tail else np.zeros((1, 0))
    Zs = []; qs = []
    for hv in itertools.product(*[range(-b[i], b[i] + 1) for i in head]):
        Z = np.zeros((len(T), n)); Z[:, tail] = T
        for i, v in zip(head, hv): Z[:, i] = v
        q = np.einsum('ij,jk,ik->i', Z, G, Z)
        k = (q <= bound * (1 + 1e-9)) & (q > 0)
        if k.any(): Zs.append(Z[k]); qs.append(q[k])
    if not Zs: return (np.zeros((0, n)), np.zeros(0)), b, size
    Z = np.vstack(Zs); q = np.concatenate(qs)
    first = np.array([row[np.nonzero(row)[0][0]] for row in Z])
    return (Z[first > 0], q[first > 0]), b, size

def check_form(A, Bc, tol=1e-6):
    """A: form, Bc: basis (columns) of the lattice.  Returns dict with own min, #min pairs, V, perfection rank."""
    n = A.shape[0]
    Lc = np.linalg.cholesky(A)
    M = Lc.T @ Bc
    U = my_lll(M)
    assert abs(int_det(U)) == 1, 'transform not unimodular'
    Uf = np.array(U, dtype=float)
    Br = Bc @ Uf
    G = Br.T @ A @ Br; G = (G + G.T) / 2
    bound = float(np.diag(G).min()) * (1 + tol)
    Z, q = fp_dfs(G, bound)
    minq = float(q.min())
    mins = np.abs(q - minq) <= tol * 10
    nmin = int(mins.sum())
    X = Z[mins] @ Br.T
    iu = np.triu_indices(n, 1)
    Vc = np.concatenate([X * X, math.sqrt(2) * X[:, iu[0]] * X[:, iu[1]]], axis=1)
    rank = int(np.linalg.matrix_rank(Vc, tol=1e-7))
    V = kappa(n) / math.sqrt(np.linalg.det(A)) / abs(np.linalg.det(Bc))
    # the lattice-free ball has radius sqrt(minq) in the A-metric: honest V uses the true minimum
    Vhonest = V * minq ** (n / 2)
    (bx, b, size) = box_enum(G, bound)
    box = None
    if bx is not None:
        Zb, qb = bx
        box = dict(minq=float(qb.min()), nmin=int((np.abs(qb - qb.min()) <= tol * 10).sum()), size=size, n_le_bound=len(qb))
    return dict(minq=minq, nmin=nmin, rank=rank, V=V, Vhonest=Vhonest, n_le_bound=len(q), box=box, boxsize=size)

if __name__ == '__main__':
    # self-test on root lattices with min 2 scaled to min 1
    def A_n(n):
        G = 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1); return G
    def D_n(n):
        G = A_n(n); G[n - 1, n - 2] = G[n - 2, n - 1] = 0; G[n - 1, n - 3] = G[n - 3, n - 1] = -1; return G
    def E_n(n):
        G = A_n(n); G[n - 1, n - 2] = G[n - 2, n - 1] = 0; G[n - 1, n - 4] = G[n - 4, n - 1] = -1; return G
    rng = np.random.default_rng(0)
    for name, G, exp_pairs in [('A2', A_n(2), 3), ('A3', A_n(3), 6), ('D4', D_n(4), 12), ('A4', A_n(4), 10), ('D5', D_n(5), 20),
                               ('E6', E_n(6), 36), ('E7', E_n(7), 63), ('E8', E_n(8), 120), ('A9', A_n(9), 45), ('D10', D_n(10), 90)]:
        n = G.shape[0]
        # disguise: random unimodular change of basis, then realize as form A on basis Bc
        Um = np.eye(n)
        for _ in range(30):
            i, j = rng.integers(0, n, 2)
            if i != j: Um[:, i] += rng.integers(-2, 3) * Um[:, j]
        G2 = Um.T @ (G / 2.0) @ Um
        Bc = np.linalg.cholesky(G2).T      # basis with Gram G2, form A = Id
        Q = np.linalg.qr(rng.standard_normal((n, n)))[0]; S = rng.standard_normal((n, n)); S = S @ S.T + n * np.eye(n)
        Ls = np.linalg.cholesky(S)          # A = S, basis = Ls^-T Q Bc  => Gram stays G2
        Bc2 = np.linalg.solve(Ls.T, Q @ Bc)
        r = check_form(S, Bc2)
        print(name, 'min', round(r['minq'], 9), 'pairs', r['nmin'], '(expected', exp_pairs, ') rank', r['rank'], 'V', round(r['V'], 4),
              'box', r['box'], 'boxsize', r['boxsize'])
