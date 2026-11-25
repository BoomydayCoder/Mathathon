"""Verifier's own tools: LLL on Gram matrices, short-vector enumeration (own C code), dual certificate."""
import numpy as np, ctypes, os, math, itertools
from scipy.optimize import nnls
_here = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(_here, 'my_enum.so'))
_lib.enum_short.restype = ctypes.c_long
_lib.enum_short.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_double, ctypes.c_long, ctypes.c_void_p,
                            ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_long)]

def kap(n):
    return math.exp((n / 2) * math.log(math.pi) - math.lgamma(n / 2 + 1))

# centre densities of the best known lattice packings (Conway-Sloane table 1.2 / Nebe-Sloane), n = 1..24
CD = {6: 1 / (8 * 3 ** 0.5), 7: 1 / 16, 8: 1 / 16, 9: 1 / (16 * 2 ** 0.5), 10: 1 / (16 * 3 ** 0.5), 11: 1 / (18 * 3 ** 0.5),
      12: 1 / 27, 13: 1 / (18 * 3 ** 0.5), 14: 1 / (16 * 3 ** 0.5), 15: 1 / (16 * 2 ** 0.5), 16: 1 / 16, 17: 1 / 16,
      18: 1 / (8 * 3 ** 0.5), 19: 1 / (8 * 2 ** 0.5), 20: 1 / 8, 21: 1 / (4 * 2 ** 0.5), 22: 1 / (2 * 3 ** 0.5), 23: 0.5, 24: 1.0,
      26: 1 / 3 ** 0.5, 28: 2.0 / 3}
def Vbest(n):
    return kap(n) * 2 ** n * CD[n]
TASK_VBEST = {6: 23.87, 7: 37.80, 8: 64.94, 9: 74.64, 10: 94.23, 11: 123.8, 12: 202.6, 13: 239.3, 14: 354.3}

def my_lll(G0, delta=0.99):
    n = G0.shape[0]
    U = np.eye(n, dtype=np.int64)
    R = np.linalg.cholesky((G0 + G0.T) / 2).T
    k = 1; it = 0
    while k < n:
        it += 1
        if it > 200000: raise RuntimeError('lll')
        for j in range(k - 1, -1, -1):
            q = int(np.rint(R[j, k] / R[j, j]))
            if q:
                U[:, k] -= q * U[:, j]
                R[:, k] -= q * R[:, j]
        if R[k, k] ** 2 + R[k - 1, k] ** 2 >= delta * R[k - 1, k - 1] ** 2:
            k += 1
        else:
            U[:, [k - 1, k]] = U[:, [k, k - 1]]
            G = U.T @ G0 @ U
            R = np.linalg.cholesky((G + G.T) / 2).T
            k = max(k - 1, 1)
    G = U.T @ G0 @ U
    return (G + G.T) / 2, U

def short(A, bound, maxc=60000, maxnodes=400_000_000):
    """all nonzero integer x (one per +-pair) with A[x] <= bound; returns X, vals (recomputed in numpy), status"""
    n = A.shape[0]
    G, U = my_lll(A)
    R = np.ascontiguousarray(np.linalg.cholesky(G).T)
    out = np.zeros((maxc, n), dtype=np.int32); vals = np.zeros(maxc)
    st = ctypes.c_int(0); nodes = ctypes.c_long(0)
    cnt = _lib.enum_short(n, R.ctypes.data, float(bound), maxc, out.ctypes.data, vals.ctypes.data, maxnodes,
                          ctypes.byref(st), ctypes.byref(nodes))
    if st.value != 0:
        raise RuntimeError('enumeration incomplete, status %d' % st.value)
    Y = out[:cnt].astype(np.int64)
    X = Y @ U.T
    v = np.einsum('ij,jk,ik->i', X, A, X) if cnt else np.zeros(0)
    return X, v, nodes.value

def min_and_pairs(A):
    X, v, _ = short(A, 1 + 1e-6)
    if len(v) == 0:
        return None, 0, X
    m = float(v.min())
    S = X[np.abs(v - 1) < 1e-8]
    return m, len(S), S

def brute_min(A, r):
    """brute force over [-r, r]^n (n small): min of A[x] over nonzero x and number of +-pairs within 1e-8 of 1"""
    n = A.shape[0]
    g = np.array(list(itertools.product(range(-r, r + 1), repeat=n)), dtype=np.int64)
    g = g[np.any(g != 0, axis=1)]
    v = np.einsum('ij,jk,ik->i', g, A, g)
    return float(v.min()), int((np.abs(v - 1) < 1e-8).sum()) // 2

def dual_certificate(Q, S):
    """lam >= 0 with sum lam_x x x^T = Q^-1 over the minimal vectors S; returns (rel residual, sum lam, #positive, min positive lam)"""
    n = Q.shape[0]
    iu = np.triu_indices(n)
    Phi = (S[:, iu[0]] * S[:, iu[1]]).astype(float)      # rows: x_i x_j (i<=j)
    c = np.linalg.inv(Q)[iu]
    lam, rn = nnls(Phi.T, c, maxiter=50 * Phi.shape[0] + 1000)
    return rn / np.linalg.norm(c), float(lam.sum()), int((lam > 1e-9).sum()), lam

def rank_perfection(S):
    n = S.shape[1]; iu = np.triu_indices(n)
    Phi = (S[:, iu[0]] * S[:, iu[1]]).astype(float)
    s = np.linalg.svd(Phi, compute_uv=False)
    return int((s > 1e-9 * s[0]).sum())
