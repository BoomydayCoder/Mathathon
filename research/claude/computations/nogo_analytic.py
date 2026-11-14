#!/usr/bin/env python3
"""nogo_analytic.py -- Theorem B (Jensen lower bound): for odd s>=3, K>=2,
 log(Delta_K(zeta(s))/K!^(s-1)) >= -E* K^2 - K*I - K log K - s*log K! + K[(s-2) log K + log(2(2pi)^(s-1)/(s-1)!) + s*e1 - e2 - 2 - log(73/16)]
with rho = sigma|[0,3/4] (density of equilibrium_s3.py), E* = int V drho - I = (9/2) log 2, I = iint log|x^2-y^2| drho drho,
e1 = int rho log x, e2 = int rho log rho.   This script: numerical values of I, int V drho, e1, e2 (mpmath; NUMERICAL),
and the bound versus exact log F(zeta(s)) (values from nogo_newton.py)."""
from mpmath import mp, mpf, quad, log, atanh, sqrt, pi, atan, factorial, loggamma
mp.dps = 25
b, c = mpf(3) / 4, mpf(5) / 4
def rho(x):
    sq = sqrt(b * b - x * x); return (3 * log((b + sq) / x) - atanh(sq / c)) / pi   # atanh(sq/b)=log((b+sq)/x)
V = lambda x: 2 * pi * x + log(1 + x * x) - 2 + 2 * x * atan(1 / x)
pts = [0, mpf(1) / 100, mpf(1) / 10, mpf(3) / 8, b]
mass = quad(rho, pts); intV = quad(lambda x: V(x) * rho(x), pts)
e1 = quad(lambda x: rho(x) * log(x), pts); e2 = quad(lambda x: rho(x) * log(rho(x)), pts)
def U(x):  # int log|x^2-y^2| rho(y) dy
    return quad(lambda y: log(abs(x * x - y * y) + mpf(10) ** -60) * rho(y), sorted(set([0, x / 2, x, (x + b) / 2, b])))
ell = 2 + 6 * log(2)
print("mass of rho          =", mass)
print("int V drho           =", intV, "  (claimed 3 log 2 - 2 =", 3 * log(2) - 2, ")")
for x in [mpf(1) / 20, mpf(1) / 5, mpf(2) / 5, mpf(3) / 5, mpf(7) / 10]:
    print(f"EL check x={float(x):.2f}: 2U(x)-V(x)+ell =", mp.nstr(2 * U(x) - V(x) + ell, 8))
mp.dps = 15
I = quad(lambda x: U(x) * rho(x), [0, mpf(1) / 50, mpf(1) / 5, mpf(1) / 2, b])
print("I = iint log|x^2-y^2| =", I, "  (claimed -2 - 1.5 log 2 =", -2 - 1.5 * log(2), ")")
Estar = intV - I
print("E* = int V - I       =", Estar, "  (9/2) log 2 =", 4.5 * log(2))
print("e1 = int rho log x   =", e1, ";  e2 = int rho log rho =", e2)
def bound(s, K):
    Ic = -2 - 1.5 * log(2)
    return (-4.5 * log(2) * K * K - K * Ic - K * log(K) - s * loggamma(K + 1)
            + K * ((s - 2) * log(K) + log(2 * (2 * pi) ** (s - 1) / factorial(s - 1)) + s * e1 - e2 - 2 - log(mpf(73) / 16)))
exact = {(7, 40): -5047.683, (7, 80): -20096.243, (7, 120): -45128.043, (5, 40): -4987.908, (9, 40): -5128.061, (3, 40): -4959.923, (11, 80): -20477.204}
print("\n s   K   Theorem-B bound   exact log F(zeta(s))   bound/K^2   exact/K^2   ok")
for (s, K), ex in sorted(exact.items()):
    bd = bound(s, K); print(f"{s:2d} {K:4d} {float(bd):14.1f} {ex:16.3f} {float(bd)/K**2:12.4f} {ex/K**2:10.4f}   {bd <= ex}")
for K in (10**3, 10**4, 10**5):
    print(f" s=7 K={K}: bound/K^2 = {float(bound(7, K))/K**2:.5f}   (-(9/2)log2 = {-4.5*float(log(2)):.5f})")
