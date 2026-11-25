"""Independent verification library for R8 (own LLL, own enumerator, own group construction)."""
import numpy as np, math, ctypes, os, sys, itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(HERE, 'myenum.so'))
_lib.my_enum.restype = ctypes.c_long
_lib.my_enum.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_double, ctypes.c_long,
                         ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_int)]

# center densities of the best known lattices (Conway-Sloane table), V = delta*kappa_n*2^n
DELTA = {2: 1 / (2 * math.sqrt(3)), 4: 1 / 8, 6: 1 / (8 * math.sqrt(3)), 8: 1 / 16, 10: 1 / (16 * math.sqrt(3)),
         12: 1 / 27, 16: 1 / 16, 20: 1 / 8, 24: 1.0}


def kappa(n):
    return math.pi ** (n / 2) / math.gamma(n / 2 + 1)


def Vbest(n):
    return DELTA[n] * kappa(n) * 2 ** n


def ldl(G):
    n = len(G)
    L = np.eye(n)
    d = np.zeros(n)
    for i in range(n):
        for j in range(i):
            L[i, j] = (G[i, j] - np.dot(L[i, :j] * d[:j], L[j, :j])) / d[j]
        d[i] = G[i, i] - np.dot(L[i, :i] ** 2, d[:i])
    return L, d


def my_lll(G0, delta=0.99):
    """G0: symmetric PD (float or int). returns unimodular int64 U with U^T G0 U LLL-reduced."""
    n = len(G0)
    Gf = np.array(G0, dtype=float)
    U = np.eye(n, dtype=np.int64)
    L = np.linalg.cholesky(Gf)
    k = 1
    guard = 0
    while k < n:
        guard += 1
        if guard > 200000:
            raise RuntimeError('LLL guard')
        for j in range(k - 1, -1, -1):
            q = int(np.rint(L[k, j] / L[j, j]))
            if q:
                U[:, k] -= q * U[:, j]
                L[k, :] -= q * L[j, :]
        if L[k, k] ** 2 + L[k, k - 1] ** 2 >= delta * L[k - 1, k - 1] ** 2:
            k += 1
        else:
            U[:, [k - 1, k]] = U[:, [k, k - 1]]
            Uf = U.astype(float)
            Gr = Uf.T @ Gf @ Uf
            L = np.linalg.cholesky((Gr + Gr.T) / 2)
            k = max(k - 1, 1)
    return U


def enum(G, bound, cap=250000, maxnodes=2_000_000_000):
    """all nonzero integer y (one per pair) with y^T G y <= bound (float G, should be reduced). returns Y (m x n int64)"""
    n = len(G)
    L, d = ldl(np.array(G, dtype=float))
    L = np.ascontiguousarray(L)
    d = np.ascontiguousarray(d)
    out = np.zeros((cap, n), dtype=np.int32)
    nodes = ctypes.c_long(0)
    st = ctypes.c_int(0)
    cnt = _lib.my_enum(n, L.ctypes.data, d.ctypes.data, float(bound), cap, out.ctypes.data, maxnodes,
                       ctypes.byref(nodes), ctypes.byref(st))
    if st.value:
        raise RuntimeError('node limit')
    if cnt > cap:
        raise RuntimeError('cap exceeded: %d' % cnt)
    return out[:cnt].astype(np.int64), nodes.value


def bareiss_det(M):
    A = [[int(x) for x in row] for row in M]
    n = len(A)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    sw = i
                    break
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def frac_inverse(M):
    n = len(M)
    A = [[Fraction(int(M[i][j])) for j in range(n)] + [Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]


def find_q(A, qmax=3000):
    """smallest q with q*A integral, A given rounded to 9 decimals"""
    for q in range(1, qmax + 1):
        Z = q * A
        if np.abs(Z - np.rint(Z)).max() < 6e-10 * q + 1e-9:
            return q
    return 0


# ---------------- groups (own construction) -----------------
def qmat(u):
    a, b, c, d = u
    return np.array([[a, -b, -c, -d], [b, a, -d, c], [c, d, a, -b], [d, -c, b, a]], dtype=float)


CYC = {5: [1, 1, 1, 1, 1], 7: [1] * 7, 8: [1, 0, 0, 0, 1], 12: [1, 0, -1, 0, 1], 13: [1] * 13}  # a_0..a_d


def companion(m):
    co = CYC[m]
    d = len(co) - 1
    M = np.zeros((d, d), dtype=np.int64)
    for j in range(d - 1):
        M[j + 1, j] = 1          # zeta * zeta^j = zeta^(j+1)
    for i in range(d):
        M[i, d - 1] = -co[i]     # zeta^d = -sum a_i zeta^i
    return M


def icosian_units():
    """120 unit icosians, exact: list of (first embedding (4,), conjugate embedding (4,)) floats"""
    s5 = math.sqrt(5)
    # values as (p,q) = p + q sqrt5
    half = (Fraction(1, 2), Fraction(0))
    one = (Fraction(1), Fraction(0))
    zero = (Fraction(0), Fraction(0))
    tau2 = (Fraction(1, 4), Fraction(1, 4))      # tau/2
    sig2 = (Fraction(1, 4), Fraction(-1, 4))     # sigma/2
    neg = lambda v: (-v[0], -v[1])
    units = set()
    for pos in range(4):
        for sgn in (1, -1):
            v = [zero] * 4
            v[pos] = one if sgn == 1 else neg(one)
            units.add(tuple(v))
    for signs in itertools.product((1, -1), repeat=4):
        units.add(tuple(half if s == 1 else neg(half) for s in signs))
    base = [zero, half, sig2, tau2]
    for perm in itertools.permutations(range(4)):
        # parity
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        if inv % 2:
            continue
        for signs in itertools.product((1, -1), repeat=3):
            vals = [zero, half if signs[0] == 1 else neg(half), sig2 if signs[1] == 1 else neg(sig2),
                    tau2 if signs[2] == 1 else neg(tau2)]
            v = [None] * 4
            for i in range(4):
                v[perm[i]] = vals[i]
            units.add(tuple(v))
    res = []
    for u in units:
        a = np.array([float(p) + float(q) * s5 for p, q in u])
        b = np.array([float(p) - float(q) * s5 for p, q in u])
        res.append((a, b))
    assert len(res) == 120, len(res)
    return res


def base_group(name):
    """returns list of ALL group elements as integer matrices on O (Z-basis coords used by eq.py), |G|, n0, extra"""
    extra = {}
    if name == 'eis':
        g = np.array([[0, -1], [1, -1]], dtype=np.int64)   # mult by omega in basis (1, omega)
        gens = [-g]
    elif name == 'gauss':
        gens = [np.array([[0, -1], [1, 0]], dtype=np.int64)]
    elif name.startswith('zeta'):
        m = int(name[4:])
        c = companion(m)
        gens = [c if m % 2 == 0 else -c]
    elif name == 'hur':
        B = np.array([[1, 0, 0, .5], [0, 1, 0, .5], [0, 0, 1, .5], [0, 0, 0, .5]])  # columns 1,i,j,(1+i+j+k)/2
        extra['covol'] = abs(np.linalg.det(B))
        us = []
        for pos in range(4):
            for s in (1, -1):
                v = [0, 0, 0, 0]
                v[pos] = s
                us.append(v)
        for signs in itertools.product((.5, -.5), repeat=4):
            us.append(list(signs))
        els = []
        for u in us:
            M = np.linalg.solve(B, qmat(u) @ B)
            assert np.abs(M - np.rint(M)).max() < 1e-9
            els.append(np.rint(M).astype(np.int64))
        return els, 24, 4, extra
    elif name == 'ico':
        sys.path.insert(0, os.path.dirname(HERE))
        import eq as their
        _, Bm = their.icosian_setup()
        extra['covol'] = abs(np.linalg.det(Bm))
        extra['Bm'] = Bm
        els = []
        coords = []
        for a, b in icosian_units():
            Lg = np.zeros((8, 8))
            Lg[:4, :4] = qmat(a)
            Lg[4:, 4:] = qmat(b)
            M = np.linalg.solve(Bm, Lg @ Bm)
            assert np.abs(M - np.rint(M)).max() < 1e-7, 'unit does not preserve their lattice'
            els.append(np.rint(M).astype(np.int64))
            c = np.linalg.solve(Bm, np.concatenate([a, b]))
            assert np.abs(c - np.rint(c)).max() < 1e-7, 'unit not in their lattice'
            coords.append(np.rint(c).astype(np.int64))
        extra['unit_coords'] = np.array(coords)
        return els, 120, 8, extra
    else:
        raise ValueError(name)
    # closure for cyclic groups
    n0 = gens[0].shape[0]
    I = np.eye(n0, dtype=np.int64)
    els = [I]
    cur = gens[0].copy()
    while not np.array_equal(cur, I):
        els.append(cur.copy())
        cur = cur @ gens[0]
        if len(els) > 1000:
            raise RuntimeError('infinite order')
    return els, len(els), n0, extra


def is_group(els):
    keys = {e.tobytes() for e in els}
    if len(keys) != len(els):
        return False
    for a in els:
        for b in els:
            if (a @ b).tobytes() not in keys:
                return False
    return True


def dimQ_character(els, k):
    """dim of G-invariants in Sym^2 of (V^k): (1/|G|) sum (chi(g)^2 + chi(g^2))/2, chi = k*tr"""
    tot = Fraction(0)
    for g in els:
        ch = k * int(np.trace(g))
        ch2 = k * int(np.trace(g @ g))
        tot += Fraction(ch * ch + ch2, 2)
    return tot / len(els)


def Qbasis_reynolds(els, k):
    n0 = els[0].shape[0]
    n = n0 * k
    gs = [np.kron(np.eye(k), g.astype(float)) for g in els]
    rows = []
    for i in range(n):
        for j in range(i, n):
            E = np.zeros((n, n))
            E[i, j] = E[j, i] = 1.0
            R = sum(g.T @ E @ g for g in gs) / len(gs)
            rows.append(R.ravel())
    R = np.array(rows)
    u, s, vt = np.linalg.svd(R, full_matrices=False)
    r = int((s > 1e-9 * s[0]).sum())
    return vt[:r].reshape(r, n, n)


def lattice_index(C):
    """C: m x n integer generators (rows) of a full-rank sublattice of Z^n; returns index (abs det of HNF)"""
    M = [[int(x) for x in row] for row in C]
    n = len(M[0])
    det = 1
    row = 0
    for c in range(n):
        while True:
            nz = [i for i in range(row, len(M)) if M[i][c] != 0]
            if len(nz) == 0:
                return 0
            p = min(nz, key=lambda i: abs(M[i][c]))
            M[row], M[p] = M[p], M[row]
            done = True
            for i in range(row + 1, len(M)):
                if M[i][c] != 0:
                    f = M[i][c] // M[row][c]
                    M[i] = [x - f * y for x, y in zip(M[i], M[row])]
                    if M[i][c] != 0:
                        done = False
            if done:
                break
        det *= abs(M[row][c])
        row += 1
    return det
