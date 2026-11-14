#!/usr/bin/env python3
"""Zudilin, 'An elementary proof of Apery's theorem' (math/0202159), Lemma 1.

R_n(t) = ((t-1)...(t-n) / (t(t+1)...(t+n)))^2,  F_n = -sum_{t>=1} R_n'(t) = u_n zeta(3) - v_n.
Exact rational arithmetic (fractions.Fraction) for all identities; mpmath (60+ digits) for F_n.
"""
from fractions import Fraction as Fr
from math import comb, lcm, factorial
import sympy as sp
import mpmath as mp

NMAX = 12
results = []


def check(name, ok, extra=""):
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))


def D(n):
    return lcm(*range(1, n + 1)) if n >= 1 else 1


def H(k, j):
    return sum((Fr(1, l ** j) for l in range(1, k + 1)), Fr(0))


def partial_fractions(n):
    """A2k = [(t+k)^2 R_n](-k), A1k = d/dt[(t+k)^2 R_n](-k), exactly, via log-derivative.
    g(t) = (t+k)^2 R_n(t) = prod_j (t-j)^2 / prod_{l!=k} (t+l)^2,
    g'/g = 2 sum_j 1/(t-j) - 2 sum_{l!=k} 1/(t+l)."""
    A2, A1 = [], []
    for k in range(n + 1):
        t = -k
        num = 1
        for j in range(1, n + 1):
            num *= (t - j)
        den = 1
        for l in range(n + 1):
            if l != k:
                den *= (t + l)
        g = Fr(num, den) ** 2
        logd = 2 * sum((Fr(1, t - j) for j in range(1, n + 1)), Fr(0)) \
            - 2 * sum((Fr(1, t + l) for l in range(n + 1) if l != k), Fr(0))
        A2.append(g)
        A1.append(g * logd)
    return A2, A1


def Rn_sym(n, t):
    num = sp.prod([t - j for j in range(1, n + 1)]) if n else sp.Integer(1)
    den = sp.prod([t + j for j in range(0, n + 1)])
    return (num / den) ** 2


print("=" * 72)
print("Lemma 1: partial fractions of R_n(t), n = 0..%d" % NMAX)
print("=" * 72)
t = sp.symbols('t')
u, v = {}, {}
for n in range(NMAX + 1):
    A2, A1 = partial_fractions(n)
    # cross-check against sympy apart for small n
    if n <= 5:
        pf = sum(sp.Rational(A2[k].numerator, A2[k].denominator) / (t + k) ** 2
                 + sp.Rational(A1[k].numerator, A1[k].denominator) / (t + k) for k in range(n + 1))
        check(f"n={n}: partial-fraction expansion equals R_n(t) (sympy cancel)",
              sp.cancel(pf - Rn_sym(n, t)) == 0)
    check(f"n={n}: A2k = C(n+k,n)^2 C(n,k)^2 for all k",
          all(A2[k] == comb(n + k, n) ** 2 * comb(n, k) ** 2 for k in range(n + 1)))
    check(f"n={n}: sum_k A1k = 0", sum(A1, Fr(0)) == 0)
    Dn = D(n)
    check(f"n={n}: D_n * A1k in Z for all k", all((Dn * a).denominator == 1 for a in A1),
          f"(lcm of denominators of A1k = {lcm(*[a.denominator for a in A1])}, D_n = {Dn})")
    # closed form A1k = 2 A2k (2H_k - H_{n+k} - H_{n-k})
    check(f"n={n}: A1k = 2 A2k (2H_k - H_(n+k) - H_(n-k))",
          all(A1[k] == 2 * A2[k] * (2 * H(k, 1) - H(n + k, 1) - H(n - k, 1)) for k in range(n + 1)))
    un = 2 * sum(A2, Fr(0))
    vn = 2 * sum((A2[k] * H(k, 3) for k in range(n + 1)), Fr(0)) + sum((A1[k] * H(k, 2) for k in range(n + 1)), Fr(0))
    u[n], v[n] = un, vn
    apery = sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))
    check(f"n={n}: u_n in Z and u_n = 2 a_n (Apery)", un.denominator == 1 and un == 2 * apery, f"u_n = {un}")
    check(f"n={n}: D_n^3 v_n in Z", (Dn ** 3 * vn).denominator == 1,
          f"v_n = {vn};  denom(v_n) = {vn.denominator}, D_n^3 = {Dn**3}, D_n^3/denom = {Fr(Dn**3, vn.denominator)}")

print()
print("A1k for n = 1, 2, 3:")
for n in (1, 2, 3):
    A2, A1 = partial_fractions(n)
    print(f"  n={n}: A2 = {[str(a) for a in A2]}, A1 = {[str(a) for a in A1]}")
print("Paper: R_1(t) = 1/t^2 + 4/(t+1)^2 - 4/t + 4/(t+1)")
A2, A1 = partial_fractions(1)
check("R_1 partial fractions match paper", A2 == [1, 4] and A1 == [-4, 4])

print()
print("Denominator table: n, denom(v_n), D_n^3, ratio D_n^3/denom(v_n), denom(2 v_n)")
for n in range(NMAX + 1):
    Dn = D(n)
    print(f"  n={n:2d}: denom(v_n) = {v[n].denominator:>14d}, D_n^3 = {Dn**3:>14d}, ratio = {Fr(Dn**3, v[n].denominator)}, denom(2v_n) = {(2*v[n]).denominator}")

print()
print("=" * 72)
print("Numerical check of F_n = -sum_{t>=1} R_n'(t) = u_n zeta(3) - v_n  (mp.dps = 70)")
print("=" * 72)
mp.mp.dps = 70


def mRprime(n, tt):
    tt = mp.mpf(tt)
    num = mp.fprod([tt - j for j in range(1, n + 1)]) if n else mp.mpf(1)
    den = mp.fprod([tt + j for j in range(0, n + 1)])
    R = (num / den) ** 2
    logd = 2 * mp.fsum([1 / (tt - j) for j in range(1, n + 1)]) - 2 * mp.fsum([1 / (tt + j) for j in range(0, n + 1)])
    return R * logd


for n in range(0, 7):
    # R_n'(t) = 0 at t = 1..n (double zeros), so sum starts at t = n+1
    Fn = -mp.nsum(lambda s: mRprime(n, s), [n + 1, mp.inf])
    pred = mp.mpf(u[n].numerator) / u[n].denominator * mp.zeta(3) - mp.mpf(v[n].numerator) / v[n].denominator
    err = abs(Fn - pred)
    check(f"n={n}: F_n numeric = u_n zeta(3) - v_n", err < mp.mpf(10) ** -55,
          f"F_n = {mp.nstr(Fn, 30)}, |diff| = {mp.nstr(err, 3)}")
check("F_0 = 2 zeta(3), F_1 = 10 zeta(3) - 12 (exact u,v)", u[0] == 2 and v[0] == 0 and u[1] == 10 and v[1] == 12)

print()
print(f"SUMMARY: {sum(results)}/{len(results)} checks PASS" + ("" if all(results) else "  -- SOME FAILED"))
