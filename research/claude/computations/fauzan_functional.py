#!/usr/bin/env python3
"""Checks of the rational functional mu_X from Fauzan's preprint "zeta(5) is irrational" (Sec. 2-3),
generalised to odd s in {3,5,7}.

f(y) = 1/(e^{2 pi y}-1),  w_s(y) = 2 y^s f^{(s-1)}(y)/(s-1)! = 2 (2pi)^{s-1} y^s sum_l l^{s-1} e^{-2 pi l y}/(s-1)!
(s=5: w = y^5 f''''/12, the paper's weight).

Exact arithmetic: fractions / sympy rationals.  Numerics: mpmath (dps 40 / 30 / 60+).
"""
from fractions import Fraction as Fr
from math import factorial, lcm, gcd, log
import random
import sympy as sp
import mpmath as mp

results = []


def check(name, ok, extra=""):
    results.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))


SLIST = (3, 5, 7)
t, x, X = sp.symbols('t x X')


def bern(k):
    """Bernoulli numbers with B_1 = -1/2 (sympy>=1.12 uses +1/2)."""
    if k == 1:
        return Fr(-1, 2)
    b = sp.bernoulli(k)
    return Fr(int(b.p), int(b.q))


def Hs(j, s):
    return sum((Fr(1, v ** s) for v in range(1, j + 1)), Fr(0))


def moment(s, e):
    """mu^{(s)}(t^e) = int_0^oo y^{2e} w_s(y) dy  (exact rational)."""
    return (-1) ** e * bern(2 * e + 2) * Fr(factorial(2 * e + s), factorial(s - 1) * factorial(2 * e + 2))


def c_s(s):
    return Fr(1, s - 1)


def pole(s, j):
    """mu^{(s)}_X(1/(t+j^2)) as (coef of X, constant)."""
    return (Fr(j ** (s - 1)), -j ** (s - 1) * Hs(j, s) + Fr(1, 2 * j) - c_s(s))


# ---------------------------------------------------------------- numerics for w_s
def eulerian(n):
    """coefficients of Eulerian polynomial A_n(q), Li_{-n}(q) = q A_n(q)/(1-q)^{n+1}."""
    A = [1]
    for m in range(1, n + 1):
        new = [0] * m
        for k in range(m):
            a_k = A[k] if k < len(A) else 0
            a_km1 = A[k - 1] if 0 <= k - 1 < len(A) else 0
            new[k] = (k + 1) * a_k + (m - k) * a_km1
        A = new
    return A


EUL = {s: eulerian(s - 1) for s in SLIST}


def w(s, y):
    y = mp.mpf(y)
    if y == 0:
        return 1 / mp.pi  # f^{(s-1)}(y) ~ (s-1)!/(2 pi y^s)
    q = mp.exp(-2 * mp.pi * y)
    omq = -mp.expm1(-2 * mp.pi * y)
    A = EUL[s]
    li = q * mp.polyval(A[::-1], q) / omq ** s  # Li_{1-s}(q)
    return 2 * (2 * mp.pi) ** (s - 1) * y ** s * li / factorial(s - 1)


def integral(s, F):
    return mp.quad(lambda y: F(y) * w(s, y), [0, 0.25, 1, 3, 6, 12, 25, mp.inf])


# ======================================================================== B1
print("=" * 76)
print("B1. Moments int_0^oo y^{2e} w_s(y) dy, s=3,5,7, e=0..6 (dps 40)")
print("=" * 76)
mp.mp.dps = 40
# sanity: series form of w_s vs closed form
ok = True
for s in SLIST:
    for y in (mp.mpf('0.3'), mp.mpf('1.7')):
        ser = 2 * (2 * mp.pi) ** (s - 1) * y ** s * mp.nsum(lambda l: l ** (s - 1) * mp.exp(-2 * mp.pi * l * y), [1, mp.inf]) / factorial(s - 1)
        ok &= abs(ser - w(s, y)) < mp.mpf(10) ** -35
check("w_s closed form (Eulerian polynomials) = series definition", ok)
print("   w_s(0+) = 1/pi for every s (paper: w -> 1/pi):", [mp.nstr(w(s, mp.mpf(10) ** -20), 10) for s in SLIST],
      mp.nstr(1 / mp.pi, 10))
for s in SLIST:
    mx = mp.mpf(0)
    mx2 = mp.mpf(0)
    for e in range(7):
        I = integral(s, lambda y: y ** (2 * e))
        f1 = 2 * mp.factorial(2 * e + s) * mp.zeta(2 * e + 2) / (mp.factorial(s - 1) * (2 * mp.pi) ** (2 * e + 2))
        q = moment(s, e)
        mx = max(mx, abs(I / f1 - 1))
        mx2 = max(mx2, abs(f1 / (mp.mpf(q.numerator) / q.denominator) - 1))
    check(f"s={s}: int y^(2e) w_s = 2(2e+s)! zeta(2e+2)/((s-1)!(2pi)^(2e+2)), e=0..6", mx < mp.mpf(10) ** -30, f"max rel err {mp.nstr(mx,3)}")
    check(f"s={s}: ... = (-1)^e B_(2e+2) (2e+s)!/((s-1)!(2e+2)!)  (exact rational)", mx2 < mp.mpf(10) ** -35, f"max rel err {mp.nstr(mx2,3)}")
    print(f"     s={s}: moments e=0..4 = {[str(moment(s, e)) for e in range(5)]}")
ok = all(moment(5, e) == (-1) ** e * bern(2 * e + 2) * Fr((2 * e + 3) * (2 * e + 4) * (2 * e + 5), 24) for e in range(30))
check("s=5 agrees with paper (2.2): (-1)^e B_{2e+2}(2e+3)(2e+4)(2e+5)/24, e=0..29", ok)

# ======================================================================== B2
print("\n" + "=" * 76)
print("B2. int_0^oo w_s(y)/(y^2+a^2) dy  vs  a^{s-1} zeta(s,a) - 1/(2a) - c_s  (dps 40)")
print("=" * 76)
mp.mp.dps = 40
for s in SLIST:
    print(f"   s={s}: residual  I - a^(s-1) zeta(s,a) + 1/(2a)  (should be -c_s):")
    resid = []
    for a in (1, 2, 3, 7, mp.mpf('2.5')):
        a = mp.mpf(a)
        I = integral(s, lambda y: 1 / (y ** 2 + a ** 2))
        r = I - a ** (s - 1) * mp.zeta(s, a) + 1 / (2 * a)
        resid.append(r)
        print(f"        a={mp.nstr(a,3):>4}: I = {mp.nstr(I, 25)}, residual = {mp.nstr(r, 25)}")
    cfit = -resid[0]
    print(f"        fitted c_{s} = {mp.nstr(cfit, 25)}  ~ identify: {mp.identify(cfit)}")
    spread = max(abs(r - resid[0]) for r in resid)
    check(f"s={s}: residual independent of a (no 1/a^k corrections), spread {mp.nstr(spread,3)}", spread < mp.mpf(10) ** -30)
    check(f"s={s}: c_s = 1/(s-1) = {c_s(s)}", abs(cfit - mp.mpf(1) / (s - 1)) < mp.mpf(10) ** -30)
    # integer a = j: j^{s-1}(zeta(s) - H_j) + 1/(2j) - c_s
    mx = mp.mpf(0)
    for j in (1, 2, 3, 7):
        I = integral(s, lambda y: 1 / (y ** 2 + j ** 2))
        cx, c0 = pole(s, j)
        val = mp.mpf(cx.numerator) / cx.denominator * mp.zeta(s) + mp.mpf(c0.numerator) / c0.denominator
        mx = max(mx, abs(I - val))
    check(f"s={s}: at a=j: I = j^(s-1)(zeta(s)-H_j^(s)) + 1/(2j) - 1/(s-1), j=1,2,3,7", mx < mp.mpf(10) ** -30)
print("   (Derivation: Hermite: zeta(s,a) = a^-s/2 + a^(1-s)/(s-1) + 2 int sin(s atan(y/a))(a^2+y^2)^(-s/2)/(e^(2pi y)-1) dy;")
print("    s-1 integrations by parts give I = 2 a^(s-1) int f(y) Im (a-iy)^(-s) dy, so c_s = 1/(s-1) exactly.)")


# ======================================================================== B3
print("\n" + "=" * 76)
print("B3. mu_X(R) direct  vs  tau_X(x^s R(-x^2)),  tau(P) = L(P^(s-2))/(s-1)!,  tau_X(1/(x-r)) = H_d(r) - X")
print("=" * 76)


def to_fr(q):
    q = sp.Rational(q)
    return Fr(int(q.p), int(q.q))


def mu_direct(s, num, den):
    """num, den: sympy expressions in t; den squarefree with roots among -j^2 (j>=1). Returns (cX, c0)."""
    N = sp.Poly(num, t, domain='QQ')
    Dp = sp.Poly(den, t, domain='QQ')
    Q, Rm = sp.div(N, Dp)
    cX, c0 = Fr(0), Fr(0)
    for (e,), co in Q.terms():
        c0 += to_fr(co) * moment(s, e)
    dD = Dp.diff(t)
    roots = sp.roots(Dp, multiple=False)
    for r, mult in roots.items():
        assert mult == 1
        j = sp.sqrt(-r)
        assert j.is_integer and j > 0, r
        res = to_fr(Rm.eval(r) / dD.eval(r))
        px, p0 = pole(s, int(j))
        cX += res * px
        c0 += res * p0
    return cX, c0


def Lfun(P):
    """L(x^k) = B_k (B_1=-1/2) on a sympy Poly in x."""
    return sum((to_fr(co) * bern(k) for (k,), co in P.terms()), Fr(0))


def tau(s, g_num, g_den):
    """tau_X on rational function g = g_num/g_den in x with simple integer poles."""
    N = sp.Poly(g_num, x, domain='QQ')
    Dp = sp.Poly(g_den, x, domain='QQ')
    Q, Rm = sp.div(N, Dp)
    c0 = Lfun(Q.diff((x, s - 2))) / factorial(s - 1) if Q.degree() >= s - 2 else Fr(0)
    cX = Fr(0)
    dD = Dp.diff(x)
    for r, mult in sp.roots(Dp, multiple=False).items():
        assert mult == 1 and r.is_integer
        res = to_fr(Rm.eval(r) / dD.eval(r))
        d = int(r) if r >= 0 else -int(r) - 1
        cX += -res
        c0 += res * Hs(d, s)
    return cX, c0


def mu_tau(s, num, den):
    gnum = sp.expand(x ** s * num.subs(t, -x ** 2))
    gden = sp.expand(den.subs(t, -x ** 2))
    return tau(s, gnum, gden)


def D(m):
    return sp.prod([t + j * j for j in range(1, m + 1)]) if m else sp.Integer(1)


tests = []
for e in range(6):
    tests.append((f"t^{e}", t ** e, sp.Integer(1)))
for j in (1, 2, 3, 5):
    tests.append((f"1/(t+{j*j})", sp.Integer(1), t + j * j))
tests.append(("t^3(t+1)^2/((t+4)(t+9)(t+25))", t ** 3 * (t + 1) ** 2, (t + 4) * (t + 9) * (t + 25)))
tests.append(("(2t^2-3t+7)/((t+1)(t+16))", 2 * t ** 2 - 3 * t + 7, (t + 1) * (t + 16)))
random.seed(1)
for _ in range(3):
    js = sorted(random.sample(range(1, 7), 3))
    nm = sum(random.randint(-5, 5) * t ** k for k in range(6))
    tests.append((f"random num deg5 / prod(t+j^2), j in {js}", nm, sp.prod([t + j * j for j in js])))
for (Nn, Kk) in ((1, 4), (2, 6)):
    for ij in (0, 2, 4, 6):
        num = D(Nn) ** 6 * t ** ij
        den = D(Kk)
        g = sp.cancel(num / den)
        nn, dd = sp.fraction(g)
        tests.append((f"D_{Nn}^6 t^{ij}/D_{Kk}", sp.expand(nn), sp.expand(dd)))

mp.mp.dps = 30
for s in SLIST:
    allok = True
    numok = True
    worst = mp.mpf(0)
    for name, num, den in tests:
        a = mu_direct(s, num, den)
        b = mu_tau(s, num, den)
        if a != b:
            allok = False
            print(f"     MISMATCH s={s} {name}: direct {a}  tau {b}")
        # numerical
        fnum = sp.lambdify(t, num, 'mpmath')
        fden = sp.lambdify(t, den, 'mpmath')
        I = integral(s, lambda y: fnum(y ** 2) / fden(y ** 2))
        with mp.workdps(120):  # X-coefficient and constant cancel heavily (up to ~25 digits)
            val = mp.mpf(a[0].numerator) / a[0].denominator * mp.zeta(s) + mp.mpf(a[1].numerator) / a[1].denominator
        err = abs(I - val) / max(1, abs(I))
        worst = max(worst, err)
        if err > mp.mpf(10) ** -20:
            numok = False
            print(f"     NUM MISMATCH s={s} {name}: integral {I}, mu {val}")
    check(f"s={s}: mu_direct == tau-pullback exactly on {len(tests)} test functions", allok)
    check(f"s={s}: mu_(zeta(s))(R) = int R(y^2) w_s(y) dy on all tests (dps 30)", numok, f"max rel err {mp.nstr(worst,3)}")
print("   example (s=5): mu_X(t^3(t+1)^2/((t+4)(t+9)(t+25))) =", mu_direct(5, t ** 3 * (t + 1) ** 2, (t + 4) * (t + 9) * (t + 25)))
print("   General-s statement: with moments (-1)^e B_{2e+2}(2e+s)!/((s-1)!(2e+2)!), pole values")
print("   j^{s-1}(X-H_j^(s)) + 1/(2j) - 1/(s-1), tau_s(P)=L(P^(s-2))/(s-1)!, tau_X(1/(x-r)) = H^(s)_d(r) - X:")
print("   mu^(s)_X(R) = tau_{s,X}(x^s R(-x^2))  (no s-dependent sign or scale).")

# ======================================================================== B4
print("\n" + "=" * 76)
print("B4. mu_(zeta(s))(R) vs -(1/(s-1)!) sum'_{x>=1} g^(s-1)(x),  g = x^s R(-x^2)")
print("    (sum' : at a pole x=r>0 of g the term is replaced by the finite part, i.e. (g - c/(x-r))^(s-1)(r))")
print("=" * 76)
mp.mp.dps = 30


def plain_sum(s, num, den):
    gx = sp.cancel(x ** s * num.subs(t, -x ** 2) / den.subs(t, -x ** 2))
    gd = sp.diff(gx, x, s - 1)
    gn, gdn = sp.fraction(sp.together(gd))
    # positive integer poles of g
    pr = [int(r) for r in sp.roots(sp.Poly(sp.fraction(gx)[1], x)).keys() if r.is_integer and r > 0]
    M = max(pr + [0]) + 1
    total = mp.mpf(0)
    for xv in range(1, M):
        if xv in pr:
            c = sp.limit((x - xv) * gx, x, xv)
            reg = sp.cancel(gx - c / (x - xv))
            fp = sp.diff(reg, x, s - 1).subs(x, xv)
            total += mp.mpf(sp.Rational(fp).p) / sp.Rational(fp).q
        else:
            v = sp.Rational(gd.subs(x, xv))
            total += mp.mpf(v.p) / v.q
    f = sp.lambdify(x, gd, 'mpmath')
    total += mp.nsum(f, [M, mp.inf])
    return -total / factorial(s - 1), gx


B4cases = []
base = (t + 1) * (t + 4) * (t + 9)
for s, es in ((3, (0, 1, 2)), (5, (0, 1, 2))):
    for e in es:
        B4cases.append((s, f"t^{e}/((t+1)(t+4)(t+9))", t ** e, base))
for s in (3, 5):
    B4cases.append((s, "1/(t+1)", sp.Integer(1), t + 1))
    B4cases.append((s, "1/((t+1)(t+4))", sp.Integer(1), (t + 1) * (t + 4)))
for s, name, num, den in B4cases:
    a = mu_direct(s, num, den)
    mu = mp.mpf(a[0].numerator) / a[0].denominator * mp.zeta(s) + mp.mpf(a[1].numerator) / a[1].denominator
    S, gx = plain_sum(s, num, den)
    gn, gdn = sp.fraction(gx)
    P = sp.div(sp.Poly(gn, x), sp.Poly(gdn, x))[0]
    degg = sp.degree(gn, x) - sp.degree(gdn, x)
    bt = Lfun(P.diff((x, s - 2))) / factorial(s - 1) if P.degree() >= s - 2 else Fr(0)
    diff = mu - S
    btm = mp.mpf(bt.numerator) / bt.denominator
    ok = abs(diff - btm) < mp.mpf(10) ** -20
    check(f"s={s}, R={name}: deg g={degg}; mu - (-(1/(s-1)!)sum') = {mp.nstr(diff, 12)} ; predicted boundary tau(poly part)={bt}", ok)
print("   Result: mu_(zeta(s))(R) = -(1/(s-1)!) sum'_{x>=1} g^(s-1)(x) + L(P^(s-2))/(s-1)!, P = polynomial part of g.")
print("   The boundary term vanishes iff deg g < s-2; for deg g = s-2 with leading coefficient c it equals c/(s-1)")
print("   (deg g >= s-1 makes the plain sum diverge).")

# ======================================================================== B5
print("\n" + "=" * 76)
print("B5. Hankel determinants Delta(X) = det[ mu^(s)_X(t^(i+j)/D_K(t)) ]_{0<=i,j<h}")
print("=" * 76)
mp.mp.dps = 400
KH = [(2, 1), (2, 2), (3, 2), (3, 3), (4, 3), (4, 4), (5, 4), (6, 5)]
rows = []
for s in SLIST:
    zs = mp.zeta(s)
    for K, h in KH:
        DK = sp.expand(D(K))
        mus = [mu_direct(s, t ** m, DK) for m in range(2 * h - 1)]
        Mat = sp.Matrix(h, h, lambda i, j: sp.Rational(mus[i + j][0].numerator, mus[i + j][0].denominator) * X +
                        sp.Rational(mus[i + j][1].numerator, mus[i + j][1].denominator))
        Delta = sp.Poly(sp.expand(Mat.det(method='berkowitz')), X)
        coeffs = [sp.Rational(c) for c in Delta.all_coeffs()]
        den_l = 1
        for c in coeffs:
            den_l = lcm(den_l, int(c.q))
        ints = [int(c * den_l) for c in coeffs]
        g = 0
        for c in ints:
            g = gcd(g, c)
        P = [c // g for c in ints]
        if P[0] < 0:
            P = [-c for c in P]
        val = mp.mpf(0)
        for c in coeffs:
            val = val * zs + mp.mpf(c.p) / c.q
        Pv = mp.mpf(0)
        for c in P:
            Pv = Pv * zs + c
        rows.append((s, K, h, Delta.degree(), val, Pv, max(abs(c) for c in P), den_l))
        # numeric Gram check for small cases
print(f" {'s':>2} {'K':>2} {'h':>2} {'deg':>3} {'Delta(zeta(s))':>14} {'log|P(zeta)|':>13} {'log max|P_i|':>13} {'log lcm den':>12}")
okpos = True
for s, K, h, dg, val, Pv, mc, dl in rows:
    okpos &= val > 0
    lp = float(mp.log(abs(Pv)))
    print(f" {s:>2} {K:>2} {h:>2} {dg:>3} {mp.nstr(val, 8):>14} {lp:>13.4f} {log(mc):>13.4f} {log(dl):>12.4f}")
check("Delta(zeta(s)) > 0 in all B5 cases (Gram positivity)", okpos)
print("   (P = primitive integer polynomial proportional to Delta, leading coeff > 0; deg = min(h, K) expected.)")

print("\n" + "=" * 76)
print(f"SUMMARY: {sum(results)}/{len(results)} checks passed")
