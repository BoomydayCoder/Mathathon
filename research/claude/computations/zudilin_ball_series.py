#!/usr/bin/env python3
"""Zudilin math/0202159: Ball's series F~_n = sum_t R~_n(t), Lemma 4 (bound), Lemma 5 (certificate (13)),
Lemma 6/7, partial fractions of R~_n.  Exact arithmetic (sympy / fractions) for identities, mpmath for numerics."""
from fractions import Fraction as Fr
from math import comb, lcm, factorial
import sympy as sp
import mpmath as mp

results = []


def check(name, ok, extra=""):
    results.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))


def D(n):
    return lcm(*range(1, n + 1)) if n >= 1 else 1


def H(k, j):
    return sum((Fr(1, l ** j) for l in range(1, k + 1)), Fr(0))


n, t = sp.symbols('n t')


def Rt_def(m, tt):
    """R~_m(t) from the product definition (8), m a concrete integer, tt sympy expr / Rational."""
    num = sp.factorial(m) ** 2 * (2 * tt + m)
    for j in range(1, m + 1):
        num *= (tt - j)
    for j in range(m + 1, 2 * m + 1):
        num *= (tt + j)
    den = sp.Integer(1)
    for j in range(0, m + 1):
        den *= (tt + j)
    return num / den ** 4


def P(nn, tt):
    return (-tt ** 6 - (8 * nn - 1) * tt ** 5 + (4 * nn ** 2 + 27 * nn + 5) * tt ** 4
            + 2 * nn * (67 * nn ** 2 + 71 * nn + 15) * tt ** 3
            + (358 * nn ** 4 + 339 * nn ** 3 + 76 * nn ** 2 - 7 * nn - 3) * tt ** 2
            + (384 * nn ** 5 + 396 * nn ** 4 + 97 * nn ** 3 - 29 * nn ** 2 - 17 * nn - 2) * tt
            + nn * (153 * nn ** 5 + 183 * nn ** 4 + 50 * nn ** 3 - 30 * nn ** 2 - 22 * nn - 4))


def St_def(m, tt):
    return Rt_def(m, tt) / ((2 * tt + m) * (tt + 2 * m - 1) * (tt + 2 * m)) * P(m, tt)


# ---------------------------------------------------------------------------
print("=" * 72)
print("(a) Reduced identity of Lemma 5 (symbolic in n, t)")
print("=" * 72)
# ratios derived by hand from (8):
r_up = (n + 1) ** 2 * (2 * t + n + 1) * (t - n - 1) * (t + 2 * n + 1) * (t + 2 * n + 2) / ((2 * t + n) * (t + n + 1) ** 5)   # R~_{n+1}/R~_n
r_down = (2 * t + n - 1) * (t + n) ** 5 / (n ** 2 * (2 * t + n) * (t - n) * (t + 2 * n - 1) * (t + 2 * n))                 # R~_{n-1}/R~_n
r_shift = (2 * t + n + 2) * t ** 5 * (t + 2 * n + 1) / ((2 * t + n) * (t - n) * (t + n + 1) ** 5)                            # R~_n(t+1)/R~_n(t)
print("  R~_(n+1)/R~_n = (n+1)^2 (2t+n+1)(t-n-1)(t+2n+1)(t+2n+2) / ((2t+n)(t+n+1)^5)")
print("  R~_(n-1)/R~_n = (2t+n-1)(t+n)^5 / (n^2 (2t+n)(t-n)(t+2n-1)(t+2n))")
print("  R~_n(t+1)/R~_n(t) = (2t+n+2) t^5 (t+2n+1) / ((2t+n)(t-n)(t+n+1)^5)")
# sanity checks of ratios against the product definition
pts = [(3, sp.Rational(7, 2)), (1, sp.Rational(5, 3)), (2, sp.Rational(-7, 11)), (5, sp.Rational(13, 2)), (4, sp.Integer(9))]
for m, tv in pts:
    ok = (Rt_def(m + 1, tv) / Rt_def(m, tv) == r_up.subs({n: m, t: tv})
          and Rt_def(m - 1, tv) / Rt_def(m, tv) == r_down.subs({n: m, t: tv})
          and Rt_def(m, tv + 1) / Rt_def(m, tv) == r_shift.subs({n: m, t: tv}))
    check(f"  ratio formulas match product definition at n={m}, t={tv} (exact)", ok)
for m in range(1, 6):
    ok = (sp.cancel(Rt_def(m + 1, t) / Rt_def(m, t) - r_up.subs(n, m)) == 0
          and sp.cancel(Rt_def(m - 1, t) / Rt_def(m, t) - r_down.subs(n, m)) == 0
          and sp.cancel(Rt_def(m, t + 1) / Rt_def(m, t) - r_shift.subs(n, m)) == 0)
    check(f"  ratio formulas match as rational functions of t, n={m}", ok)

Sred = P(n, t) / ((2 * t + n) * (t + 2 * n - 1) * (t + 2 * n))                       # S~_n(t)/R~_n(t)
Sred_shift = r_shift * P(n, t + 1) / ((2 * t + n + 2) * (t + 2 * n) * (t + 2 * n + 1))  # S~_n(t+1)/R~_n(t)
lhs = (n + 1) ** 3 * r_up - (2 * n + 1) * (17 * n ** 2 + 17 * n + 5) + n ** 3 * r_down
rhs = Sred_shift - Sred
diff = sp.cancel(sp.together(lhs - rhs))
check("reduced identity (14)/R~_n holds identically in (n,t)", diff == 0, "" if diff == 0 else f"residual numerator = {sp.factor(sp.numer(sp.together(diff)))}")
for m in range(1, 7):
    d = ((m + 1) ** 3 * Rt_def(m + 1, t) - (2 * m + 1) * (17 * m ** 2 + 17 * m + 5) * Rt_def(m, t) + m ** 3 * Rt_def(m - 1, t)
         - (St_def(m, t + 1) - St_def(m, t)))
    check(f"  identity (14) exact as rational function of t, n={m}", sp.cancel(sp.together(d)) == 0)

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("(b) S~_n(1) = 0 for n>=1, S~_n(t) = O(t^-2), R~_n(t) = O(t^-5)")
print("=" * 72)
for m in range(1, 9):
    S = sp.cancel(St_def(m, t))
    s1 = S.subs(t, 1)
    degS = sp.degree(sp.numer(S), t) - sp.degree(sp.denom(S), t)
    Rr = sp.cancel(Rt_def(m, t))
    degR = sp.degree(sp.numer(Rr), t) - sp.degree(sp.denom(Rr), t)
    check(f"n={m}: S~_n(1) = 0, ord_inf S~_n = {degS} (<= -2), ord_inf R~_n = {degR} (<= -5)",
          s1 == 0 and degS <= -2 and degR <= -5)
S0 = sp.cancel(St_def(0, t))
print(f"  info: S~_0(t) = {sp.factor(S0)}, S~_0(1) = {S0.subs(t, 1)}")

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("(e)/(c) Exact partial fractions R~_n = sum_{j=1..4} sum_k C_jk/(t+k)^j")
print("=" * 72)


def ser_mul(a, b, N=4):
    c = [Fr(0)] * N
    for i in range(N):
        if a[i] == 0:
            continue
        for j in range(N - i):
            c[i + j] += a[i] * b[j]
    return c


def pf_tilde(m, N=4):
    """Laurent coefficients at t=-k: (t+k)^4 R~_m(t) = g0 + g1 e + g2 e^2 + g3 e^3 + ..., e = t+k."""
    C = {j: [] for j in range(1, 5)}
    for k in range(m + 1):
        g = [Fr(factorial(m) ** 2)] + [Fr(0)] * (N - 1)
        lin = [(m - 2 * k, 2)] + [(-k - j, 1) for j in range(1, m + 1)] + [(-k + j, 1) for j in range(m + 1, 2 * m + 1)]
        for c0, c1 in lin:
            g = ser_mul(g, [Fr(c0), Fr(c1)] + [Fr(0)] * (N - 2), N)
        for l in range(m + 1):
            if l == k:
                continue
            c = Fr(l - k)
            inv = [Fr((-1) ** i * comb(i + 3, 3)) / c ** (4 + i) for i in range(N)]
            g = ser_mul(g, inv, N)
        C[4].append(g[0]); C[3].append(g[1]); C[2].append(g[2]); C[1].append(g[3])
    return C


def uv_apery(m):
    A2 = [comb(m + k, m) ** 2 * comb(m, k) ** 2 for k in range(m + 1)]
    A1 = [2 * A2[k] * (2 * H(k, 1) - H(m + k, 1) - H(m - k, 1)) for k in range(m + 1)]
    return Fr(2 * sum(A2)), 2 * sum((A2[k] * H(k, 3) for k in range(m + 1)), Fr(0)) + sum((A1[k] * H(k, 2) for k in range(m + 1)), Fr(0))


for m in (0, 1):
    C = pf_tilde(m)
    expr = sum(sp.Rational(C[j][k].numerator, C[j][k].denominator) / (t + k) ** j for j in range(1, 5) for k in range(m + 1))
    print(f"  R~_{m}(t) = {expr}")
C1 = pf_tilde(1)
paper1 = {4: [-2, 2], 3: [5, 5], 2: [-5, 5], 1: [0, 0]}
check("R~_1 = -2/t^4 + 2/(t+1)^4 + 5/t^3 + 5/(t+1)^3 - 5/t^2 + 5/(t+1)^2 (paper's Lemma 7 display)",
      all(C1[j] == paper1[j] for j in range(1, 5)))
check("R~_0 = 2/t^3", pf_tilde(0)[3] == [2] and pf_tilde(0)[4] == [0] and pf_tilde(0)[2] == [0] and pf_tilde(0)[1] == [0])

ut, vt = {}, {}
for m in range(0, 9):
    C = pf_tilde(m)
    if m <= 4:
        expr = sum(sp.Rational(C[j][k].numerator, C[j][k].denominator) / (t + k) ** j for j in range(1, 5) for k in range(m + 1))
        check(f"n={m}: partial fractions reproduce R~_n (sympy cancel)", sp.cancel(expr - Rt_def(m, t)) == 0)
    s = {j: sum(C[j], Fr(0)) for j in range(1, 5)}
    u_ap, v_ap = uv_apery(m)
    vv = sum((C[j][k] * H(k, j) for j in range(1, 5) for k in range(m + 1)), Fr(0))
    ut[m], vt[m] = s[3], vv
    check(f"n={m}: sum_k C1k = 0, zeta(2)-coeff sum_k C2k = 0, zeta(4)-coeff sum_k C4k = 0", s[1] == 0 and s[2] == 0 and s[4] == 0,
          f"(sums: C1={s[1]}, C2={s[2]}, C4={s[4]})")
    check(f"n={m}: u~_n = sum_k C3k equals u_n = 2 a_n", s[3] == u_ap, f"u~_n = {s[3]}")
    check(f"n={m}: v~_n = sum_jk C_jk H_k^(j) equals v_n exactly (=> F~_n = F_n exactly)", vv == v_ap)
    Dn = D(m)
    dC = [lcm(*[C[j][k].denominator for k in range(m + 1)]) for j in range(1, 5)]
    print(f"     denominators: lcm denom(C_1k..C_4k) = {dC};  denom(u~_n) = {s[3].denominator};  denom(v~_n) = {vv.denominator};"
          f"  D_n = {Dn}, D_n^3 = {Dn**3}, D_n^4 = {Dn**4};  D_n^3 v~_n in Z: {(Dn**3*vv).denominator == 1};  D_n^4 v~_n in Z: {(Dn**4*vv).denominator == 1}")
    # individual-term integrality as in Ball's argument: D_n^{4-j} C_jk in Z ?
    okj = all((Dn ** (4 - j) * C[j][k]).denominator == 1 for j in range(1, 5) for k in range(m + 1))
    print(f"     D_n^(4-j) C_jk in Z for all j,k: {okj}")

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("(c)/(d) Numerics: F~_n = sum_{t>=n+1} R~_n(t) vs u_n zeta(3) - v_n   (mp.dps = 90)")
print("=" * 72)
mp.mp.dps = 90


def mRt(m, x):
    x = mp.mpf(x)
    num = mp.mpf(factorial(m)) ** 2 * (2 * x + m)
    for j in range(1, m + 1):
        num *= (x - j)
    for j in range(m + 1, 2 * m + 1):
        num *= (x + j)
    den = mp.fprod([x + j for j in range(0, m + 1)])
    return num / den ** 4


z3 = mp.zeta(3)
q = (mp.sqrt(2) - 1) ** 4
Ft = {}
for m in range(0, 9):
    Ft[m] = mp.nsum(lambda x: mRt(m, x), [m + 1, mp.inf])
    pred = mp.mpf(ut[m].numerator) / ut[m].denominator * z3 - mp.mpf(vt[m].numerator) / vt[m].denominator
    err = abs(Ft[m] - pred)
    rel = err / abs(Ft[m])
    bound = 20 * (m + 1) ** 4 * q ** m
    ratio = Ft[m] / Ft[m - 1] if m >= 1 else None
    check(f"n={m}: F~_n = u_n zeta(3) - v_n (rel err {mp.nstr(rel, 3)})", rel < mp.mpf(10) ** -50,
          f"F~_n = {mp.nstr(Ft[m], 25)}" + (f", F~_n/F~_(n-1) = {mp.nstr(ratio, 12)}" if ratio else ""))
    check(f"n={m}: 0 < F~_n < 20(n+1)^4 (sqrt2-1)^(4n)", 0 < Ft[m] < bound, f"(bound = {mp.nstr(bound, 8)}, F~/bound = {mp.nstr(Ft[m]/bound, 6)})")
print(f"  (sqrt2-1)^4 = {mp.nstr(q, 20)}")
check("F~_0 = 2 zeta(3) and F~_1 = 10 zeta(3) - 12 numerically", abs(Ft[0] - 2 * z3) < mp.mpf(10) ** -80 and abs(Ft[1] - (10 * z3 - 12)) < mp.mpf(10) ** -80)

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("(f) Lemma 4 analysis")
print("=" * 72)
tau = sp.symbols('tau')
x = tau + sp.Rational(1, 2)
lhs = sp.expand(tau ** 5 * (tau + 2) - (tau - 1) * (tau + 1) ** 5)
print(f"  tau^5(tau+2) - (tau-1)(tau+1)^5 = {lhs}   (degree {sp.degree(lhs, tau)}, not 6)")
cand_user = -(2 * x ** 4 - 5 * x ** 2 - sp.Rational(7, 8))
cand_paper = -x * (2 * x ** 4 - 5 * x ** 2 - sp.Rational(7, 8))
check("user's candidate  -[2(tau+1/2)^4 - 5(tau+1/2)^2 - 7/8]  equals it", sp.expand(lhs - cand_user) == 0)
check("corrected identity: tau^5(tau+2)-(tau-1)(tau+1)^5 = -(tau+1/2)[2(tau+1/2)^4 - 5(tau+1/2)^2 - 7/8]",
      sp.expand(lhs - cand_paper) == 0)
print(f"  factorization: {sp.factor(lhs)}")
tau0 = -sp.Rational(1, 2) + sp.sqrt(sp.Rational(5, 4) + sp.sqrt(2))
check("tau0 = -1/2 + sqrt(5/4 + sqrt2) is a root (exact)", sp.simplify(sp.expand(lhs.subs(tau, tau0))) == 0)
realroots = [r for r in sp.Poly(lhs, tau).nroots(n=30) if abs(sp.im(r)) < 1e-25]
print(f"  real roots: {[sp.N(sp.re(r), 20) for r in realroots]};  tau0 = {sp.N(tau0, 25)}")
check("tau0 is the unique real root > 1", sum(1 for r in realroots if sp.re(r) > 1) == 1 and sp.N(tau0) > 1)

mp.mp.dps = 60
T0 = -mp.mpf(1) / 2 + mp.sqrt(mp.mpf(5) / 4 + mp.sqrt(2))
f = lambda T: 5 * T * mp.log(T) + (T + 2) * mp.log(T + 2) - (T - 1) * mp.log(T - 1) - 5 * (T + 1) * mp.log(T + 1)
fp = lambda T: mp.log(T ** 5 * (T + 2) / ((T - 1) * (T + 1) ** 5))
target = 4 * mp.log(mp.sqrt(2) - 1)
check("f'(tau0) = 0", abs(fp(T0)) < mp.mpf(10) ** -55, f"f'(tau0) = {mp.nstr(fp(T0), 5)}")
check("f'(tau) = numeric derivative of f at tau=2.7", abs(mp.diff(f, mp.mpf('2.7')) - fp(mp.mpf('2.7'))) < mp.mpf(10) ** -40)
alt = 2 * mp.log(T0 + 2) + mp.log(T0 - 1) - 5 * mp.log(T0 + 1)
check("f(tau0) = 2log(tau0+2) + log(tau0-1) - 5log(tau0+1)", abs(f(T0) - alt) < mp.mpf(10) ** -55)
check("f(tau0) = 4 log(sqrt2 - 1)", abs(f(T0) - target) < mp.mpf(10) ** -55, f"f(tau0) = {mp.nstr(f(T0), 30)}, 4log(sqrt2-1) = {mp.nstr(target, 30)}")
grid = [1 + mp.mpf(i) / 100 for i in range(1, 2000)]
check("sup_{tau>1} f = f(tau0) (grid tau in (1,21))", max(f(g) for g in grid) <= f(T0) + mp.mpf(10) ** -50)

C = mp.e ** 2 * (2 * mp.zeta(5) + 5 * mp.zeta(4) + 2 * mp.zeta(3))
print(f"  e^2 (2zeta5 + 5zeta4 + 2zeta3) = {mp.nstr(C, 20)}")
check("literal constant claim e^2(2zeta(5)+5zeta(4)+2zeta(3)) < 20", C < 20)
# what is really needed: e^2 (n+1)^2 (2zeta5 + 5n zeta4 + 2n^2 zeta3) < 20 (n+1)^4, i.e. coefficientwise in n
c2, c1, c0 = 20 - 2 * mp.e ** 2 * mp.zeta(3), 40 - 5 * mp.e ** 2 * mp.zeta(4), 20 - 2 * mp.e ** 2 * mp.zeta(5)
print(f"  20(n+1)^2 - e^2(2zeta5 + 5n zeta4 + 2n^2 zeta3) = {mp.nstr(c2, 8)} n^2 + {mp.nstr(c1, 8)} n + {mp.nstr(c0, 8)}")
check("all coefficients positive => e^2(2zeta5+5n zeta4+2n^2 zeta3) < 20(n+1)^2 for all n>=0", c2 > 0 and c1 > 0 and c0 > 0)
d2, d0 = 20 - 2 * mp.e ** 2 * mp.zeta(5), 20 - 2 * mp.e ** 2 * mp.zeta(3)
check("swapped version e^2(2zeta3+5n zeta4+2n^2 zeta5) < 20(n+1)^2 for all n>=0 (coefficients positive)", d2 > 0 and c1 > 0 and d0 > 0)
# the termwise bound used: (2t+n)(t+2n)/(t^4(t+n)) <= 2/t^3 + 5n/t^4 + 2n^2/t^5 (natural) and <= 2/t^5+5n/t^4+2n^2/t^3 (as printed) for t>=n+1, n>=1
okA = okB = True
for m in range(1, 31):
    for tt in range(m + 1, 20 * m + 1):
        T = Fr(tt)
        val = Fr((2 * tt + m) * (tt + 2 * m), tt ** 4 * (tt + m))
        okA &= val < 2 / T ** 3 + 5 * m / T ** 4 + 2 * m * m / T ** 5
        okB &= val <= 2 / T ** 5 + 5 * m / T ** 4 + 2 * m * m / T ** 3
check("(2t+n)(t+2n)/(t^4(t+n)) < 2/t^3+5n/t^4+2n^2/t^5, n=1..30, t=n+1..20n (exact)", okA)
check("(2t+n)(t+2n)/(t^4(t+n)) <= 2/t^5+5n/t^4+2n^2/t^3 (paper's printed order), n=1..30, t=n+1..20n (exact)", okB)
print(f"  27(sqrt2-1)^4 = {mp.nstr(27 * (mp.sqrt(2) - 1) ** 4, 20)}")
check("27(sqrt2-1)^4 = 0.7948... < 1", mp.nstr(27 * (mp.sqrt(2) - 1) ** 4, 4) == '0.7948' and 27 * (mp.sqrt(2) - 1) ** 4 < 1)

# inequality (12)
mp.mp.dps = 40
worst = mp.mpf(0); wpt = None; ok12 = True
for m in range(1, 31):
    rhs = mp.e ** 2 * (m + 1) ** 2 * q ** m
    for tt in range(m + 1, 20 * m + 1):
        lhsv = mRt(m, tt) * mp.mpf(tt) ** 4 * (tt + m) / ((2 * tt + m) * (tt + 2 * m))
        r = lhsv / rhs
        if r > worst:
            worst, wpt = r, (m, tt)
        if not lhsv < rhs:
            ok12 = False
check("inequality (12): R~_n(t) t^4(t+n)/((2t+n)(t+2n)) < e^2(n+1)^2(sqrt2-1)^(4n), n=1..30, t=n+1..20n", ok12,
      f"(max LHS/RHS = {mp.nstr(worst, 8)} at (n,t) = {wpt})")
# the sharper intermediate bound of (11): LHS/(n+1)^2 < (1+1/n)^{2n} e^{n f(t/n)}
ok11 = True
for m in range(1, 31):
    for tt in range(m + 1, 20 * m + 1):
        lhsv = mRt(m, tt) * mp.mpf(tt) ** 4 * (tt + m) / ((2 * tt + m) * (tt + 2 * m) * (m + 1) ** 2)
        ff = f(mp.mpf(tt) / m)
        if not lhsv < (1 + mp.mpf(1) / m) ** (2 * m) * mp.e ** (m * ff):
            ok11 = False
check("intermediate bound (11): LHS < (1+1/n)^(2n) e^(n f(t/n)), n=1..30, t=n+1..20n", ok11)

print()
print(f"SUMMARY: {sum(results)}/{len(results)} checks PASS" + ("" if all(results) else "  -- SOME FAILED (see above; some FAILs are expected findings)"))
