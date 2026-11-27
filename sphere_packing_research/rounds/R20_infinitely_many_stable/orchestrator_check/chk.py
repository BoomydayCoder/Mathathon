# Independent check of the R20 ideal model: uses only build() from the agent's code (lattice construction),
# then my own exact tests and my own (pure Python, exact rational) enumeration of short vectors.
import os, sys, time
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ideal_model as im

def bareiss(M):
    A = [list(map(int, r)) for r in M]; n = len(A); prev = 1; sign = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None: return 0
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]

def lll_gram(G, delta=Fraction(99, 100)):
    """exact LLL on an integer Gram matrix; returns (reduced Gram, transform T with G' = T G T^t)."""
    n = len(G)
    T = [[int(i == j) for j in range(n)] for i in range(n)]
    G = [list(map(int, r)) for r in G]
    def gso():
        mu = [[Fraction(0)] * n for _ in range(n)]; B = [Fraction(0)] * n
        for i in range(n):
            for j in range(i + 1):
                s = Fraction(G[i][j])
                for l in range(j):
                    s -= mu[i][l] * mu[j][l] * B[l]
                if j < i: mu[i][j] = s / B[j]
                else: B[i] = s
        return mu, B
    def addrow(i, j, q):   # b_i -= q b_j
        for c in range(n): T[i][c] -= q * T[j][c]
        for c in range(n): G[i][c] -= q * G[j][c]
        for c in range(n): G[c][i] -= q * G[c][j]
    k = 1
    mu, B = gso()
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                addrow(k, j, q); mu, B = gso()
        if B[k] >= (delta - mu[k][k - 1] ** 2) * B[k - 1]:
            k += 1
        else:
            T[k], T[k - 1] = T[k - 1], T[k]
            G[k], G[k - 1] = G[k - 1], G[k]
            for r in G: r[k], r[k - 1] = r[k - 1], r[k]
            mu, B = gso(); k = max(k - 1, 1)
    return G, T

def enum_short(G, bound):
    """all x != 0 (up to sign) with x G x^t <= bound; floating GSO with small slack, exact recheck."""
    n = len(G)
    mu = [[0.0] * n for _ in range(n)]; B = [0.0] * n
    for i in range(n):
        for j in range(i + 1):
            s = float(G[i][j])
            for l in range(j): s -= mu[i][l] * mu[j][l] * B[l]
            if j < i: mu[i][j] = s / B[j]
            else: B[i] = s
    sols = []; x = [0] * n; import math
    def rec(i, rem):
        if i < 0:
            if any(x):
                sols.append(x[:])
            return
        c = -sum(x[j] * mu[j][i] for j in range(i + 1, n))
        rad = math.sqrt(max(rem, 0.0) / B[i]) + 1e-9
        lo = math.ceil(c - rad); hi = math.floor(c + rad)
        for v in range(lo, hi + 1):
            x[i] = v
            d = (v - c) ** 2 * B[i]
            if d <= rem + 1e-7: rec(i - 1, rem - d)
        x[i] = 0
    rec(n - 1, bound + 1e-6)
    out = []
    for v in sols:
        # keep one of +-v
        fz = next(t for t in v if t)
        if fz > 0: out.append(v)
    return out

def qf(M, v):
    n = len(v); return sum(v[i] * int(M[i][j]) * v[j] for i in range(n) if v[i] for j in range(n) if v[j])

# usage: chk.py [p ...]   (default: 13 37; the file out_p61_p277_p661.txt was produced with 61 277 661)
PRIMES = [int(a) for a in sys.argv[1:]] or [13, 37]
for (p, N) in [(p, 12) for p in PRIMES]:
    t0 = time.time()
    G, E, Bz, F = im.build(p, N, verbose=False)
    n = len(G)
    G = [[int(x) for x in r] for r in G]; E = [[int(x) for x in r] for r in E]
    assert all(G[i][j] == G[j][i] for i in range(n) for j in range(n))
    even = all(G[i][i] % 2 == 0 for i in range(n))
    dG = bareiss(G); dE = bareiss(E)
    # eta = G^{-1} E must satisfy eta^2 = p  <=>  E G^{-1} E = p G ; check exactly via adjugate-free solve: E * X = p*G with X = G^{-1}E
    # exact: solve G X = E over Q
    from fractions import Fraction as Fr
    A = [[Fr(G[i][j]) for j in range(n)] + [Fr(E[i][j]) for j in range(n)] for i in range(n)]
    for c in range(n):
        pv = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[pv] = A[pv], A[c]
        inv = 1 / A[c][c]; A[c] = [a * inv for a in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    X = [row[n:] for row in A]                       # eta in lattice coordinates (acting on columns)
    X2 = [[sum(X[i][l] * X[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
    eta_sq_ok = all(X2[i][j] == (p if i == j else 0) for i in range(n) for j in range(n))
    eta_int = all(x.denominator == 1 for r in X for x in r)
    print('p=%d N=%d rank=%d even=%s detG=%d detE=p^%s eta integral=%s eta^2=p:%s  [%.0fs]' % (
        p, N, n, even, dG, (len(str(abs(dE))) and __import__('math').log(abs(dE), p)), eta_int, eta_sq_ok, time.time() - t0), flush=True)
    Gr, T = lll_gram(G)
    print('   LLL diag:', sorted(Gr[i][i] for i in range(n)), '[%.0fs]' % (time.time() - t0), flush=True)
    roots = enum_short(Gr, 2)
    # back to original coordinates: v_orig = v_red * T
    from collections import Counter
    cnt = Counter()
    for v in roots:
        w = [sum(v[i] * T[i][j] for i in range(n)) for j in range(n)]
        nrm = qf(G, w); s2 = qf(E, w)      # s = w.E.w / 2
        assert nrm == 2 and s2 % 2 == 0
        s = abs(s2 // 2)
        Nnu4 = p * 1 - s * s               # 4*N(nu') = p*l^2 - s^2 with norm 2l = 2
        cnt[(s, Fraction(Nnu4, 4))] += 1
    print('   roots (pairs): %d ; (|s|, N(nu\')) -> pairs: %s ; threshold (N/4)^2 = %d  [%.0fs]' % (
        len(roots), dict(cnt), (N // 4) ** 2, time.time() - t0), flush=True)
