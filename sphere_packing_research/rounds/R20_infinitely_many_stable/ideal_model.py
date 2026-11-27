"""R20: explicit 'Craig-type' ideal-lattice model of the SL_2(p) Weil lattice W_p (p = 5 mod 8) and of its
split-torus fixed sublattices (quotient-curve lattices), in EXACT integer arithmetic.

Model (see proof.md, Theorem 1):  N | p-1 with k=(p-1)/N odd, N = 4m' (m' odd).  M_N = subfield of Q(zeta_p) of degree N,
L = M_N(i) (degree 2N), P1,P2 the two primes of L over p (i = iota mod P1), lam = N_{Q(zeta_p)/M_N}(1-zeta) (generator of p_N).
   L_N(p) := ( P1^{(m'+1)/2} P2^{-(m'-1)/2} ,  Tr_{L/Q}(x xbar)/(2p) )          (rank 2N, even unimodular)
For N = p-1 this is (claimed) the even unimodular SL_2(p)-lattice W_p; for N < p-1 it is Fix(C_k)/sqrt(k).
Coordinates: z = x*lam^{(m'-1)/2} in P1^{m'} subset O_{M_N}[i],  z = y1 + i y2,  y = sum a_j theta_j (Gaussian periods),
   z in P1^{m'}  <=>  sum_j (a1_j + iota a2_j) gamma^{j l} = 0 mod p,  l = 0..m'-1   (gamma = g^k of order N)
   form = (1/p) Tr_{M_N/Q}((y1 y1bar' + y2 y2bar') u^{-(m'-1)/2}),  u = lam*lambar.
Also E = Gram of (x, sqrt(p) y); s(w) = w^T E w / 2;  N(nu(w)) = (p (w.w)^2 - 4 s^2)/16 ; stability threshold m'^2.
"""
import sys, math, time, json, collections
from fractions import Fraction
import numpy as np
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', 'R2_voronoi_descent'))
import voronoi as V
# R20: the round-2 enumerator has a fixed coordinate buffer int x[64]; for rank > 64 use a local copy with x[1024]
import ctypes as _ct, os as _os
_big = _ct.CDLL(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'fp_enum_big.so'))
_big.fp_enum.restype = _ct.c_long
_big.fp_enum.argtypes = V._lib.fp_enum.argtypes
V._lib = _big


def primitive_root(p):
    fac = set(); n = p - 1; d = 2
    while d * d <= n:
        while n % d == 0:
            fac.add(d); n //= d
        d += 1
    if n > 1: fac.add(n)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g


def fund_unit(p):
    """fundamental unit (a + b sqrt p)/2 of Q(sqrt p), p = 1 mod 4, by the continued fraction of (1+sqrt p)/2 (PQa);
    returns (a, b, sgn) with a^2 - p b^2 = sgn in {-4, +4}."""
    P, Q = 1, 2
    A2, A1 = 0, 1; B2, B1 = 1, 0; G2, G1 = -1, 2
    sq = math.isqrt(p); i = 0
    while True:
        a = (P + sq) // Q
        A2, A1 = A1, a * A1 + A2; B2, B1 = B1, a * B1 + B2; G2, G1 = G1, a * G1 + G2
        P = a * Q - P; Q = (p - P * P) // Q
        if Q == 2:
            sgn = G1 * G1 - p * B1 * B1
            assert sgn in (-4, 4), sgn
            return G1, B1, sgn
        i += 1


class Field:
    """M_N subset Q(zeta_p), basis = Gaussian periods theta_j = sum_{c in g^j H} zeta^c, j = 0..N-1."""
    def __init__(self, p, N):
        assert (p - 1) % N == 0
        self.p, self.N = p, N; self.k = k = (p - 1) // N
        g = self.g = primitive_root(p)
        self.cos = [0] * p                       # coset index of c != 0
        e = 1
        for j in range(p - 1):
            self.cos[e] = j % N; e = e * g % p
        self.cosets = [[c for c in range(1, p) if self.cos[c] == j] for j in range(N)]
        # theta_i * theta_j = k*[C_j = -C_i] * 1 + sum_l n(i,j;l) theta_l,   1 = -sum theta_l
        M = np.zeros((N, N, N), dtype=np.int64)
        reps = [pow(g, l, p) for l in range(N)]
        for i in range(N):
            for l in range(N):
                gl = reps[l]
                for c in self.cosets[i]:
                    d = (gl - c) % p
                    if d: M[i, self.cos[d], l] += 1
        half = (N // 2) % N
        for i in range(N):
            M[i, (i + half) % N, :] -= k          # k * 1 = -k * sum theta_l
        self.M = M
        self.conj = [(j + half) % N for j in range(N)]

    def mul(self, a, b):
        N = self.N; out = [0] * N; M = self.M
        for i in range(N):
            if a[i] == 0: continue
            for j in range(N):
                if b[j] == 0: continue
                ab = a[i] * b[j]
                row = M[i, j]
                for l in range(N):
                    if row[l]: out[l] += ab * int(row[l])
        return out

    def from_poly(self, c):
        """c = coefficient list (length p) of an element of Z[X]/(X^p-1) that is constant on cosets."""
        N = self.N
        for j in range(N):
            assert len({c[x] for x in self.cosets[j]}) == 1
        return [c[self.cosets[j][0]] - c[0] for j in range(N)]

    def u_elem(self):
        """u = prod_{h in H} (1-zeta^h)(1-zeta^-h) as period coordinates."""
        p = self.p
        c = [0] * p; c[0] = 1
        for h in self.cosets[0]:
            c = [2 * c[x] - c[(x - h) % p] - c[(x + h) % p] for x in range(p)]
        return self.from_poly(c)

    def inverse(self, a):
        N = self.N
        # matrix of multiplication by a: column j = a*theta_j
        A = [[Fraction(0)] * (N + 1) for _ in range(N)]
        for j in range(N):
            ej = [0] * N; ej[j] = 1
            col = self.mul(a, ej)
            for i in range(N): A[i][j] = Fraction(col[i])
        for i in range(N): A[i][N] = Fraction(-1)     # 1 = -sum theta
        for c in range(N):
            piv = next(r for r in range(c, N) if A[r][c] != 0)
            A[c], A[piv] = A[piv], A[c]
            inv = 1 / A[c][c]
            A[c] = [x * inv for x in A[c]]
            for r in range(N):
                if r != c and A[r][c] != 0:
                    f = A[r][c]
                    A[r] = [x - f * y for x, y in zip(A[r], A[c])]
        return [A[i][N] for i in range(N)]

    def sqrtp(self):
        return [(-1) ** j for j in range(self.N)]

    def trace(self, a):
        return -sum(a)


def build(p, N, verbose=True, variant='L1'):
    """integer Gram matrices (G, E) of L_N(p) and the basis matrix."""
    t0 = time.time()
    assert p % 8 == 5 and (p - 1) % N == 0 and ((p - 1) // N) % 2 == 1 and N % 8 == 4
    F = Field(p, N); k = F.k; mp = N // 4; e2 = (mp - 1) // 2
    u = F.u_elem()
    one = [-1] * N
    c = [Fraction(x) for x in one]
    if variant == 'L0':        # L_0-type: ideal P1^{-(3m'-1)/2} P2^{-(m'-1)/2}, form Tr(eps x xbar/(2 sqrt p))
        e2 = (3 * mp - 1) // 2
    if e2:
        ui = F.inverse(u)
        for _ in range(e2): c = F.mul(c, ui)
    if variant == 'L0':
        a_, b_, sg_ = fund_unit(p)
        assert sg_ == -4
        epsp = [Fraction(-b_ * p, 2) + Fraction(a_, 2) * sq for sq in F.sqrtp()]     # eps*sqrt p = (b p + a sqrt p)/2, 1 = -sum theta
        c = F.mul(c, epsp)
    cs = F.mul(c, F.sqrtp())
    # functionals phi(y) = Tr(y c), phiE(y) = Tr(y c sqrt p)
    phi = []; phiE = []
    for l in range(N):
        el = [0] * N; el[l] = 1
        phi.append(F.trace(F.mul(c, el))); phiE.append(F.trace(F.mul(cs, el)))
    P = [[None] * N for _ in range(N)]; PE = [[None] * N for _ in range(N)]
    for j in range(N):
        for l in range(N):
            row = F.M[j, F.conj[l]]
            P[j][l] = sum(int(row[t]) * phi[t] for t in range(N)) / p
            PE[j][l] = sum(int(row[t]) * phiE[t] for t in range(N)) / p
    # construction A: sum_j (a1_j + iota a2_j) gamma^{jl} = 0 mod p, l < mp
    gam = pow(F.g, k, p)
    iota = pow(F.g, (p - 1) // 4, p)
    H = [[pow(gam, j * l, p) for j in range(N)] + [iota * pow(gam, j * l, p) % p for j in range(N)] for l in range(mp)]
    n = 2 * N
    # row reduce H mod p
    Hm = [r[:] for r in H]; piv = []
    r = 0
    for col in range(n):
        pr = next((i for i in range(r, mp) if Hm[i][col] % p), None)
        if pr is None: continue
        Hm[r], Hm[pr] = Hm[pr], Hm[r]
        inv = pow(Hm[r][col], p - 2, p)
        Hm[r] = [x * inv % p for x in Hm[r]]
        for i in range(mp):
            if i != r and Hm[i][col] % p:
                f = Hm[i][col]
                Hm[i] = [(x - f * y) % p for x, y in zip(Hm[i], Hm[r])]
        piv.append(col); r += 1
        if r == mp: break
    assert r == mp
    free = [cidx for cidx in range(n) if cidx not in piv]
    Bz = []
    for f in free:
        v = [0] * n; v[f] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-Hm[i][f]) % p
            if v[pc] > p // 2: v[pc] -= p
        Bz.append(v)
    for pc in piv:
        v = [0] * n; v[pc] = p; Bz.append(v)
    # Gram
    def gram(Pm):
        G = [[None] * n for _ in range(n)]
        # block diag(P,P)
        BP = []
        for v in Bz:
            a1 = v[:N]; a2 = v[N:]
            w1 = [sum(a1[j] * Pm[j][l] for j in range(N) if a1[j]) for l in range(N)]
            w2 = [sum(a2[j] * Pm[j][l] for j in range(N) if a2[j]) for l in range(N)]
            BP.append(w1 + w2)
        for i in range(n):
            for j in range(n):
                G[i][j] = sum(BP[i][t] * Bz[j][t] for t in range(n) if Bz[j][t])
        return G
    G = gram(P); E = gram(PE)
    assert all(x.denominator == 1 for row in G for x in row), 'G not integral'
    assert all(x.denominator == 1 for row in E for x in row), 'E not integral'
    Gi = np.array([[int(x) for x in row] for row in G], dtype=object)
    Ei = np.array([[int(x) for x in row] for row in E], dtype=object)
    assert (Gi == Gi.T).all() and (Ei == Ei.T).all()
    even = all(int(Gi[i, i]) % 2 == 0 for i in range(n))
    if verbose:
        print('p=%d N=%d k=%d m\'=%d rank=%d: integral Gram, even=%s  [%.1fs]' % (p, N, k, mp, n, even, time.time() - t0), flush=True)
    return Gi, Ei, np.array(Bz, dtype=object), F


def exact_det(G):
    """Bareiss determinant of an integer (object) matrix."""
    A = [list(map(int, r)) for r in G]; n = len(A); prev = 1; sign = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None: return 0
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[-1][-1]


def analyze(p, N, tstep=0.1, do_circle=True, maxnorm=4, out=None):
    Gi, Ei, Bz, F = build(p, N)
    n = 2 * N; mp = N // 4
    det = exact_det(Gi) if n <= 80 else None
    Gf = np.array(Gi, dtype=float); Ef = np.array(Ei, dtype=float)
    # check E = G * eta with eta^2 = p:  (G^-1 E)^2 = p
    Mx = np.linalg.solve(Gf, Ef)
    eta_ok = float(np.abs(Mx @ Mx - p * np.eye(n)).max())
    print('   det = %s, |(G^-1E)^2 - p| = %.1e' % (det, eta_ok), flush=True)
    res = {'p': p, 'N': N, 'rank': n, 'det': str(det), 'even': bool(all(int(Gi[i, i]) % 2 == 0 for i in range(n)))}
    thr = mp * mp
    t1 = time.time()
    Xs, vals, flag = V.short_vectors(Gf, maxnorm + 0.01, maxc=400000, maxnodes=2_000_000_000)
    Xo = np.array(Xs, dtype=object)
    cnt = collections.Counter(); minN = None; viol = 0
    for x in Xo:
        nrm = int(x @ Gi @ x); s2 = int(x @ Ei @ x)        # s = s2/2
        assert s2 % 2 == 0
        s = s2 // 2
        num = p * nrm * nrm - 4 * s * s
        assert num % 16 == 0, (nrm, s)
        Nnu = num // 16
        cnt[(nrm, abs(s), Nnu)] += 1
        if minN is None or Nnu < minN: minN = Nnu
        if Nnu < thr: viol += 1
    desc = ', '.join('norm %d |s|=%d N(nu)=%d: %d pairs' % (kk[0], kk[1], kk[2], v) for kk, v in sorted(cnt.items()))
    print('   vectors of norm <= %d (flag %d): %s' % (maxnorm, flag, desc[:1500]))
    print('   threshold m\'^2 = %d ; min N(nu) among them = %s ; violators = %d   [%.1fs]' % (thr, minN, viol, time.time() - t1), flush=True)
    res.update({'short': {str(kk): v for kk, v in cnt.items()}, 'thr': thr, 'minN_short': minN, 'viol_short': viol, 'flag': int(flag)})
    if do_circle:
        a, b, sgn = fund_unit(p); eps = (a + b * math.sqrt(p)) / 2; L = math.log(eps)
        Gp = (Gf + Ef / math.sqrt(p)) / 2; Gm = (Gf - Ef / math.sqrt(p)) / 2
        theta = 4.0 * mp * mp / p
        bound = 2 * math.sqrt(theta) * math.cosh(tstep / 2) * (1 + 1e-9)
        best = None; fail = 0; tot = 0; t1 = time.time(); seen = {}
        for t in np.arange(-L, L + 1e-9, tstep):
            Gt = math.exp(t) * Gp + math.exp(-t) * Gm
            try:
                Xs, vals, flag = V.short_vectors(Gt, bound, maxc=300000, maxnodes=500_000_000)
            except Exception as ex:
                fail += 1; continue
            if flag == 2: fail += 1
            for x in np.array(Xs, dtype=object):
                tot += 1
                nrm = int(x @ Gi @ x); s = int(x @ Ei @ x) // 2
                Nnu = (p * nrm * nrm - 4 * s * s) // 16
                seen[Nnu] = seen.get(Nnu, 0) + 1
                if best is None or Nnu < best[0]: best = (Nnu, nrm, s, float(t))
        print('   t-circle (log eps = %.3f, %s): %d vectors below bound, min N(nu) = %s (threshold %d), N-values seen: %s  [%.1fs]' % (
            L, 'exhaustive' if fail == 0 else 'INCOMPLETE %d failures' % fail, tot, best, thr, dict(sorted(seen.items())[:8]), time.time() - t1), flush=True)
        res.update({'circle_best': best, 'circle_fail': fail, 'log_eps': L})
    if out is not None: out.append(res)
    return res


if __name__ == '__main__':
    out = []
    for arg in sys.argv[1:]:
        p, N = map(int, arg.split(':'))
        try:
            analyze(p, N, out=out)
        except Exception as ex:
            print('p=%d N=%d FAILED: %r' % (p, N, ex), flush=True)
    json.dump(out, open(os.path.join(_HERE, 'ideal_model_results_%s.json') % '_'.join(a.replace(':', 'x') for a in sys.argv[1:4]), 'w'), indent=1, default=str)
