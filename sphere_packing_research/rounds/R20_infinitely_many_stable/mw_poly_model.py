"""R20: elementary 'polynomial' form of the Mordell-Weil lattice (proof.md, Corollary 4.2), checked numerically.
p = 5 mod 8, m=(p-1)/4, iota^2=-1 mod p.  Lattice of classes (c,c') in Z^p x Z^p modulo (1,0),(0,1) with
  c = f(k), c' = f'(k) mod p,  deg f, f' <= (3m-1)/2,  deg(f - iota f') <= (m-1)/2;
  V = |c-mean|^2+|c'-mean|^2,  S = c^T Q c + c'^T Q c'  (Q = Paley matrix ((j-l)/p)),
  deg = V/(2p),  N(q_F) = (V^2 - S^2/p)/(16p^2),  threshold ((p-1)/4)^2."""
import sys, math, collections
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, '../R10_weil_lattices/orchestrator')
import ideal_model as IM, proto as P
V_ = IM.V

def hnf_basis(rows):
    """integer row lattice basis (rows: list of int lists), simple Euclidean elimination."""
    A = [r[:] for r in rows]; n = len(A[0]); basis = []; col = 0
    for col in range(n):
        piv = [r for r in A if r[col] != 0]
        rest = [r for r in A if r[col] == 0]
        while len(piv) > 1:
            piv.sort(key=lambda r: abs(r[col]))
            b = piv[0]
            new = [b]
            for r in piv[1:]:
                q = r[col] // b[col]
                r2 = [x - q * y for x, y in zip(r, b)]
                if r2[col] != 0: new.append(r2)
                elif any(r2): rest.append(r2)
            piv = new
        if piv: basis.append(piv[0])
        A = rest
    return basis

def run(p, sign=1, variant=0):
    m = (p - 1) // 4; g = IM.primitive_root(p); iota = pow(g, m, p)
    if sign < 0: iota = p - iota
    d1 = (3 * m - 1) // 2; d2 = (m - 1) // 2
    gens = []
    for e in range(d1 + 1):                       # f' = k^e, f = iota f'
        gens.append([iota * pow(k, e, p) % p for k in range(p)] + [pow(k, e, p) for k in range(p)])
    for e in range(d2 + 1):                       # f = k^e, f' = 0  (so f - iota f' has degree <= d2)
        gens.append([pow(k, e, p) for k in range(p)] + [0] * p)
    for j in range(2 * p):
        v = [0] * (2 * p); v[j] = p; gens.append(v)
    # difference coordinates (kills the two all-ones directions)
    def diff(v):
        return [v[k] - v[0] for k in range(1, p)] + [v[p + k] - v[p] for k in range(1, p)]
    B = hnf_basis([diff(v) for v in gens])
    n = len(B); assert n == 2 * (p - 1), n
    # Gram of V (centered square sum) and S (Paley form) on representatives (0, d_1.., 0, d'_1..)
    leg = [0] + [1 if pow(a, (p - 1) // 2, p) == 1 else -1 for a in range(1, p)]
    Q = np.array([[leg[(j - l) % p] for l in range(p)] for j in range(p)], dtype=object)
    C = np.eye(p, dtype=object) * p - np.ones((p, p), dtype=object)          # p * centering
    rep = np.array([[0] + b[:p - 1] + [0] + b[p - 1:] for b in B], dtype=object)
    r1 = rep[:, :p]; r2 = rep[:, p:]
    GV = (r1 @ C @ r1.T + r2 @ C @ r2.T)            # = p * V-form
    GS = (r1 @ Q @ r1.T + r2 @ Q @ r2.T)
    # 2 deg = V/p = GV/p^2 : check integrality and determinant p^(p-1)
    assert all(int(x) % (p * p) == 0 for x in GV.flat), 'V/p not integral'
    G2 = np.array([[int(x) // (p * p) for x in row] for row in GV], dtype=object)      # Gram of 2*deg
    Gr, U = P.lll_gram_fraction([[int(x) for x in row] for row in G2]); U = np.array(U, dtype=object).T
    G2r = U.T @ G2 @ U; GSr = U.T @ GS @ U
    det = IM.exact_det(G2r)
    k = 0; d = det
    while d % p == 0: d //= p; k += 1
    print('p=%d: rank %d, Gram(2deg) integral, det = p^%d * %d (expected p^%d), even diag: %s' % (p, n, k, d, p - 1, all(int(G2r[i, i]) % 2 == 0 for i in range(n))), flush=True)
    bound = 2 * ((p - 1) // 2) + 2
    Xs, vals, flag = V_.short_vectors(np.array(G2r, dtype=float), bound + 0.01, maxc=600000, maxnodes=2_000_000_000)
    cnt = collections.Counter()
    for x in np.array(Xs, dtype=object):
        twodeg = int(x @ G2r @ x); S = int(x @ GSr @ x)
        Vv = twodeg * p                      # V = 2 deg * p
        num = Vv * Vv * p - S * S            # p*(V^2 - S^2/p)
        assert num % (16 * p ** 3) == 0, (twodeg, S)
        cnt[(twodeg // 2, num // (16 * p ** 3))] += 1
    print('   maps of degree <= %d (pairs +-): ' % (bound // 2) + ', '.join('deg %d N(q_F)=%d: %d' % (a, b, c) for (a, b), c in sorted(cnt.items())) + '   [threshold N >= %d, genus %d]' % (m * m, (p - 1) // 2), flush=True)

if __name__ == '__main__':
    for a in sys.argv[1:]:
        run(int(a))
