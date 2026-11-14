#!/usr/bin/env python3
"""
hankel_referee_checks.py -- independent referee checks (s = 3) for
  (I)  the analytic bound  ||p_n||^2 <= C n^2 theta^{4n} n!^4 / K!^2,  p_n(t) = prod_{j<=n}(t-(theta j)^2),
       via the one-variable inequality (H_{theta,delta}) on ALL of r > 0;
  (II) the small-prime bound L8 / (B2) and the referee's variant
          v_p(c_{i,i'}) >= -2 sum_k N_{p^k}(i,i') - 6 floor(log_p 2K) - [p=2].
Standalone (python-flint, numpy).  Prints PASS/FAIL per check.  Does not import the draft's code except
in check (I.4), which compares against the draft's exact norms (hankel_hadamard.norms_exact).

Referee's proof of (I) in brief (see final report for details):
  ratio(y) := p_n(y^2)^2 w(y)/D_K(y^2) / (theta^{4n} n!^4/K!^2) = F_n(y/theta)^2 w(y)/E_K(y),
  F_n(z) = prod_{j<=n}|1-z^2/j^2|,  E_K(y) = prod_{j<=K}(1+y^2/j^2).
  (a) log F_n(z) <= n phi(z/n) + 1 + log 2 + log n    (monotonicity of u -> log|1-z^2/u^2| on (0,z) and (z,inf))
  (b) log E_K(y) >= J_K(y) - log(1+y^2) - 2 >= J_n(y) - log(1+y^2) - 2,  J_c(y) = int_0^c log(1+y^2/u^2) du
  (c) w(y) <= (1+2 pi y)^3 e^{-2 pi y} / pi
  (H) on r = z/n = y/(theta n) > 0  gives  2 n phi(r) - 2 pi y - J_n(y) <= -delta y.
"""
import sys, math, time
import numpy as np
from flint import arb, ctx, fmpq, fmpq_poly, fmpz_poly, fmpq_mat

RESULTS = []
def report(name, ok, info=""):
    RESULTS.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {info}")
    sys.stdout.flush()

# ----------------------------------------------------------------------------------------------
# (I.1) certify (H_{theta,delta}) on (0, infinity)
# ----------------------------------------------------------------------------------------------
def H_certify(theta_q, delta_q, eps_q=fmpq(1, 1000), R=60, tol_depth=40):
    ctx.prec = 100
    th = arb(theta_q.p) / theta_q.q; de = arb(delta_q.p) / delta_q.q; pi = arb.pi()
    # (0,1]: phi(r) <= 2 log2 r^2 (positive power series, phi(r)/r^2 increasing);
    # RHS >= 2 pi th r + 2 th r arctan(1/th)  (log term >= 0; arctan(1/(th r)) >= arctan(1/th) for r <= 1)
    # => RHS - LHS >= r [th(2 pi + 2 arctan(1/th) - de) - 4 log2 r] >= r [th(2pi+2arctan(1/th)-de) - 4 log 2]
    c01 = th * (2 * pi + 2 * (1 / th).atan() - de) - 4 * arb(2).log()
    ok01 = bool(c01 > 0)
    # [1, 1+eps]: (r+1)log(r+1) <= (2+eps)log(2+eps); -(r-1)log(r-1) <= eps log(1/eps) (eps<1/e);
    # every RHS term is nondecreasing in r (x arctan(1/x) is increasing), so RHS >= RHS(1).
    e = arb(eps_q.p) / eps_q.q
    lhs = 2 * ((2 + e) * (2 + e).log() + e * (1 / e).log()) + de * th * (1 + e)
    rhs = 2 * pi * th + (1 + th * th).log() + 2 * th * (1 / th).atan()
    c1 = rhs - lhs
    ok1 = bool(c1 > 0)
    # [1+eps, R]: adaptive bisection with ball arithmetic
    def g(lo, hi):
        r = arb((lo + hi) / 2, (hi - lo) / 2)
        phi = (r + 1) * (r + 1).log() - (r - 1) * (r - 1).log()
        return 2 * pi * th * r + (1 + th * th * r * r).log() + 2 * th * r * (1 / (th * r)).atan() - 2 * phi - de * th * r
    lo0 = fmpq(1) + eps_q
    stack = [(lo0, fmpq(R), 0)]
    okmid = True; nboxes = 0; minlow = None
    while stack:
        lo, hi, d = stack.pop()
        v = g(arb(lo.p) / lo.q, arb(hi.p) / hi.q)
        if v > 0:
            nboxes += 1
            low = float(v.mid()) - float(v.rad())
            if minlow is None or low < minlow[0]:
                minlow = (low, float(lo))
            continue
        if d > tol_depth:
            okmid = False; break
        m = (lo + hi) / 2
        stack.append((lo, m, d + 1)); stack.append((m, hi, d + 1))
    # r >= R: phi(r) = 2 + int_{r-1}^{r+1} log x dx <= 2 + 2 log(r+1); RHS >= 2 pi th r;
    # need (2pi - de) th r - 4 log(r+1) - 4 >= 0 at r = R and nondecreasing for r >= R.
    Rb = arb(R)
    okR = bool((2 * pi - de) * th * Rb - 4 * (Rb + 1).log() - 4 > 0) and bool((2 * pi - de) * th - 4 / (Rb + 1) > 0)
    return ok01, ok1, okmid, okR, c01, c1, nboxes, minlow

def H_minslack_float(theta, delta):
    r = np.linspace(1.0 + 1e-12, 60, 3_000_001)
    phi = (r + 1) * np.log(r + 1) - (r - 1) * np.log(r - 1)
    d = 2 * np.pi * theta * r + np.log1p(theta ** 2 * r ** 2) + 2 * theta * r * np.arctan(1 / (theta * r)) - 2 * phi - delta * theta * r
    k = np.argmin(d)
    return d[k], r[k]

# ----------------------------------------------------------------------------------------------
# (I.2) numerics for the integrand ratio
# ----------------------------------------------------------------------------------------------
def logw(y):
    a = np.pi * y
    return np.log(4 * np.pi ** 2) + 3 * np.log(y) - 2 * a + np.log1p(np.exp(-2 * a)) - 3 * np.log(-np.expm1(-2 * a))

def ratio_scan(K, theta, delta=0.05):
    h = theta / 300
    y = np.arange(h / 2, theta * 3 * K + 60, h)
    base = logw(y)
    for j in range(1, K + 1):
        base -= np.log1p((y / j) ** 2)
    acc = np.zeros_like(y); out = []
    for n in range(K):
        if n >= 1:
            with np.errstate(divide='ignore'):
                acc += 2 * np.log(np.abs(1 - (y / (theta * n)) ** 2))
        lr = acc + base
        sup = lr.max()
        integ = math.log(np.sum(np.exp(lr - sup)) * h) + sup
        bound = 4 + 2 * math.log(2) - math.log(math.pi) + 2 * math.log(max(n, 1)) + 3 * np.log1p(2 * np.pi * y) + np.log1p(y * y) - delta * y
        out.append((n, sup, integ, float(np.max(lr - bound))))
    return out

# ----------------------------------------------------------------------------------------------
# (II) exact entries of the two-sided Newton matrix, computed independently
# ----------------------------------------------------------------------------------------------
def vp_int(n, p):
    n = int(n)
    if n == 0: return 10 ** 9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v
def vp(q, p):
    q = fmpq(q)
    return 10 ** 9 if q == 0 else vp_int(q.p, p) - vp_int(q.q, p)
def flog(x, p):  # floor(log_p x) for integer x >= 1
    k, q = 0, p
    while q <= x: k += 1; q *= p
    return k
def primes_upto(n):
    return [q for q in range(2, n + 1) if all(q % r for r in range(2, int(q ** 0.5) + 1))]

def mu_mon(e):
    return fmpq((-1) ** e) * fmpq.bernoulli(2 * e + 2) * fmpq(2 * e + 3, 2)

def tau(Q):  # Q: fmpq_poly in x; tau(Q) = L(Q')/2, L(x^k) = B_k (B_1 = -1/2)
    cs = Q.derivative().coeffs()
    return sum((cs[k] * fmpq.bernoulli(k) for k in range(len(cs))), fmpq(0)) / 2

def entries(K):
    H = [fmpq(0)]
    for v in range(1, K + 1): H.append(H[-1] + fmpq(1, v ** 3))
    m_a = [None] + [-fmpq(j * j) * H[j] - fmpq(1, 2) + fmpq(1, 2 * j) for j in range(1, K + 1)]
    m_b = [None] + [fmpq(j * j) for j in range(1, K + 1)]
    mus = [mu_mon(e) for e in range(K + 1)]
    C = {}
    for i in range(K):
        for ip in range(K):
            A, B = i, K - ip
            if A < B:  # pole entry: sum_j rho_j m_j
                a = fmpq(0); b = fmpq(0)
                for j in range(A + 1, B + 1):
                    den = 1
                    for m in range(A + 1, B + 1):
                        if m != j: den *= (m * m - j * j)
                    a += m_a[j] / den; b += m_b[j] / den
                C[(i, ip)] = ('pole', a, b, A, B)
            else:      # polynomial entry: mu(prod_{j=B+1}^{A}(t+j^2))
                P = fmpz_poly([1])
                for j in range(B + 1, A + 1): P *= fmpz_poly([j * j, 1])
                cs = P.coeffs()
                C[(i, ip)] = ('poly', sum((fmpq(int(cs[k])) * mus[k] for k in range(len(cs))), fmpq(0)), fmpq(0), A, B)
    return C

def check_II(K, do_det):
    C = entries(K)
    ok_B2 = ok_var = ok_pole_int = ok_poly_int = True
    minslack_B2 = minslack_var = 10 ** 9
    if do_det:
        # independent determinant: det(Al + X Be) and the monomial Hankel matrix via partial fractions
        pass
    for p in primes_upto(3 * K + 5):
        L3 = flog(3 * K, p); L2 = flog(2 * K, p); L1 = flog(K, p)
        for (i, ip), (kind, a, b, A, B) in C.items():
            v = min(vp(a, p), vp(b, p))
            if v >= 10 ** 8: continue
            SN = 0; q = p
            while q <= 3 * K:
                SN += (B - 1) // q - (A + 1) // q; q *= p
            b2 = -7 * L3 - 2 * SN
            var = -2 * SN - 6 * L2 - (1 if p == 2 else 0)
            if v < b2: ok_B2 = False
            if v < var: ok_var = False
            minslack_B2 = min(minslack_B2, v - b2); minslack_var = min(minslack_var, v - var)
            if kind == 'pole':
                if v < -2 * SN - 5 * L1 - L2 - (1 if p == 2 else 0): ok_pole_int = False
            else:
                d = A - B
                vd = sum(d // p ** k for k in range(1, 40))
                if v < 2 * vd - 2 * L2: ok_poly_int = False
    return ok_B2, ok_var, ok_pole_int, ok_poly_int, minslack_B2, minslack_var

def det_check(K):
    """det of Newton matrix == det of monomial Hankel matrix G_K, both as polys in X (evaluated at X = 0..K)."""
    C = entries(K)
    H = [fmpq(0)]
    for v in range(1, K + 1): H.append(H[-1] + fmpq(1, v ** 3))
    D = fmpz_poly([1])
    for j in range(1, K + 1): D *= fmpz_poly([j * j, 1])
    Dp = D.derivative()
    mus = [mu_mon(e) for e in range(2 * K)]
    def mom(e, X):  # mu_X(t^e / D_K)
        tot = fmpq(0)
        num = fmpz_poly([0] * e + [1])
        if e >= K:
            qq, _ = divmod(num, D); cs = qq.coeffs()
            tot += sum((fmpq(int(cs[k])) * mus[k] for k in range(len(cs))), fmpq(0))
        for j in range(1, K + 1):
            rj = fmpq((-j * j) ** e, int(Dp(-j * j)))
            tot += rj * (fmpq(j * j) * (X - H[j]) - fmpq(1, 2) + fmpq(1, 2 * j))
        return tot
    vals = []
    ok = True
    pvals = {}
    for X in range(0, K + 1):
        G = fmpq_mat(K, K, [mom(i + j, fmpq(X)) for i in range(K) for j in range(K)])
        N = fmpq_mat(K, K, [C[(i, ip)][1] + X * C[(i, ip)][2] for i in range(K) for ip in range(K)])
        if G.det() != N.det(): ok = False
        vals.append(N.det())
    # interpolate to get coefficients, then Gauss valuations vs the bounds
    xs = list(range(K + 1)); dd = vals[:]
    n = len(xs)
    for k in range(1, n):
        for i in range(n - 1, k - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) / fmpq(xs[i] - xs[i - k])
    poly = fmpq_poly([0])
    for k in range(n - 1, -1, -1):
        poly = poly * fmpq_poly([-xs[k], 1]) + fmpq_poly([dd[k]])
    cs = [c for c in poly.coeffs() if c != 0]
    okb = True; rows = []
    for p in primes_upto(3 * K + 5):
        v = min(vp(c, p) for c in cs)
        b = -K * (6 * flog(2 * K, p) + (1 if p == 2 else 0))
        if v < b: okb = False
        if p <= 7: rows.append((p, v, b, -7 * K * flog(3 * K, p)))
    return ok, okb, rows

# ----------------------------------------------------------------------------------------------
if __name__ == '__main__':
    t0 = time.time()
    print("=== (I.1) certification of (H_{theta,delta}) on (0,inf) ===")
    for (t, d, expect) in [(fmpq(2, 5), fmpq(1, 20), True), (fmpq(21, 50), fmpq(2, 5), True), (fmpq(9, 20), fmpq(9, 10), True)]:
        ok01, ok1, okmid, okR, c01, c1, nb, ml = H_certify(t, d)
        ms, rs = H_minslack_float(float(t), float(d))
        report(f"(H) theta={t} delta={d}", ok01 and ok1 and okmid and okR,
               f"(0,1]: const {float(c01.mid()):.4f}>0 {ok01}; [1,1.001]: {float(c1.mid()):.4f}>0 {ok1}; [1.001,60]: {okmid} ({nb} boxes, min certified lower bound {ml[0]:.5f} near r={ml[1]:.3f}); r>=60: {okR}; float min slack {ms:.6f} at r={rs:.4f}")
    for th in (0.38, 0.39):
        ms, rs = H_minslack_float(th, 0.0)
        report(f"(H) FAILS for theta={th}, delta=0 (expected negative)", ms < 0, f"min RHS-LHS on r>=1: {ms:.5f} at r={rs:.4f}")
    print("=== (I.2) bound (c) for w ===")
    y = np.linspace(1e-6, 30, 300001)
    lhs = logw(y); rhs = 3 * np.log1p(2 * np.pi * y) - 2 * np.pi * y - math.log(math.pi)
    report("w(y) <= (1+2 pi y)^3 e^{-2 pi y}/pi on grid", bool(np.all(lhs <= rhs + 1e-12)), f"max(log lhs - log rhs) = {np.max(lhs - rhs):.3e}")
    print("=== (I.3) integrand ratio p_n^2 omega_K / (theta^{4n} n!^4/K!^2): max over n<K ===")
    for theta in (0.40, 0.42, 0.45, 0.38):
        print(f"  theta={theta}:   K | max_n log sup_y ratio | max_n log int ratio | /log K | max_n,y (log ratio - referee bound)")
        for K in (5, 10, 20, 40, 80, 160):
            o = ratio_scan(K, theta)
            ms = max(x[1] for x in o); mi = max(x[2] for x in o); mv = max(x[3] for x in o)
            print(f"        {K:4d} | {ms:8.3f} | {mi:8.3f} | {mi / math.log(K):7.3f} | {mv:8.3f}")
            if theta >= 0.4:
                RESULTS.append(mv <= 0)
    print("=== (I.4) quadrature vs exact norms (draft's exact moments), K=20, theta=2/5 ===")
    try:
        sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
        from hankel_hadamard import norms_exact
        K = 20; nr, _, _ = norms_exact(K, 2, 5); o = ratio_scan(K, 0.4)
        ctx.prec = 6000 + 80 * K * K; z = arb(3).zeta(); md = 0
        for n, (a, b) in enumerate(nr):
            v = arb(a.p) / arb(a.q) + z * arb(b.p) / arb(b.q)
            tgt = (arb(2) / 5) ** (4 * n) * arb.fac_ui(n) ** 4 / arb.fac_ui(K) ** 2
            md = max(md, abs(float((v / tgt).log().mid()) - o[n][2]))
        ctx.prec = 100
        report("quadrature == exact mu_{zeta(3)}(p_n^2/D_K)", md < 1e-8, f"max log-difference {md:.2e}")
    except Exception as ex:
        report("quadrature vs exact", False, repr(ex))
    print("=== (II.1) tau(binom(x,k)) = (-1)^{k-1} H_k/(k+1), k <= 40 ===")
    ok = True; Bx = fmpq_poly([1]); Hk = fmpq(0)
    for k in range(1, 41):
        Bx = Bx * fmpq_poly([fmpq(-(k - 1), k), fmpq(1, k)]); Hk += fmpq(1, k)
        if tau(Bx) != fmpq((-1) ** (k - 1)) * Hk / (k + 1): ok = False
    report("tau(binom(x,k)) formula, k=1..40", ok)
    ok = all(tau(fmpq_poly([0, 0, 0] + [0] * (2 * e) + [1]) * 1) * (-1) ** e == mu_mon(e) for e in range(30))
    report("mu(t^e) = tau(x^3 (-x^2)^e), e<30", ok)
    print("=== (II.2) entrywise bounds, all primes p <= 3K+5 ===")
    for K in (8, 12, 20, 30, 45, 64):
        okB2, okv, okp, okq, s2, sv = check_II(K, False)
        report(f"K={K}: (B2) / variant / pole-intermediate / poly-intermediate", okB2 and okv and okp and okq,
               f"{okB2} {okv} {okp} {okq}; min slack (B2) {s2}, variant {sv}")
    print("=== (II.3) det(Newton) == det(G_K) and v_p(det) >= -K(6 floor(log_p 2K) + [p=2]) ===")
    for K in (8, 12, 16):
        ok, okb, rows = det_check(K)
        report(f"K={K}: det equality, determinant bound", ok and okb, f"{ok} {okb}; (p, v_p(det), referee bound, draft L8 bound): {rows}")
    print(f"time {time.time() - t0:.1f}s")
    print("ALL PASS" if all(RESULTS) else "SOME FAIL")
