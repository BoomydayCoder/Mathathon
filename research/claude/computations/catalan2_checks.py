#!/usr/bin/env python3
"""
catalan2_checks.py -- NUMERICAL checks (mpmath, 30 digits) of the formulas used in catalan2_lib.py / the report.

 A. pole values at mixed nodes:  int_0^inf w_B(y)/(y^2+a^2) dy = a zeta(2,a) - 1/(2a) - 1,
                                 int_0^inf w_I(y)/(y^2+a^2) dy = a eta(2,a) - 1/(2a),
    zeta(2,1/4) = pi^2+8G, zeta(2,3/4) = pi^2-8G, eta(2,1) = pi^2/12, eta(2,1/2) = 4G; monomial moments.
 B. the pi^2-free test function for w_B on S = {1/4,3/4,...} (K nodes): Psi = sum_a r_a/(t+a^2), r_a = 1/(a D'(-a^2)):
    B1  mu_B(t^e Psi) computed by quadrature = exact c0 + cG*G (cP = 0) for e <= K-2; for e = K-1 the pi^2 part is back;
    B2  Psi(t) = (1/pi) int_0^inf s^{-1/2} [1/D(t) - 1/D(s)]/(s-t) ds   (> 0: 1/D is decreasing);
    B3  lemma: sum_a r_a a = 0  =>  int_0^inf (y^2 Psi(y^2) - c) dy = 0 with c = sum_a r_a = lim t Psi(t)  (so c > 0).
"""
import sys
import mpmath as mp
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from catalan2_lib import *

mp.mp.dps = 30
pi, G = mp.pi, mp.catalan
Q = lambda x: mp.mpf(int(x.p)) / mp.mpf(int(x.q))
BR = [0, 0.25, 0.5, 1, 2, 4, 8, 16, 40]


def wB(y):
    e = mp.exp(2 * pi * y)
    return 4 * pi * y ** 2 * e / (e - 1) ** 2


def wI(y):
    return pi * y ** 2 * mp.cosh(pi * y) / mp.sinh(pi * y) ** 2


worst = mp.mpf(0)
for kernel, w, nodes in (('bose', wB, [fmpq(k, 4) for k in range(1, 13)]), ('sinh', wI, [fmpq(k, 2) for k in range(1, 9)])):
    for a in nodes:
        c0, cG, cP = pole_value(kernel, a)
        lhs = mp.quad(lambda y: w(y) / (y * y + Q(a) ** 2), BR)
        worst = max(worst, abs(lhs - (Q(c0) + Q(cG) * G + Q(cP) * pi ** 2)))
    for e in range(4):
        worst = max(worst, abs(mp.quad(lambda y: w(y) * y ** (2 * e), BR) - Q(mu_mono(kernel, e))))
print("A. max |quadrature - closed form| over 20 pole values and 8 monomial moments:", mp.nstr(worst, 3))

K = 8
S = [fmpq(2 * k + 1, 4) for k in range(K)]
dp = dprime(S)
res = {a: 1 / (a * dp[a]) for a in S}
mom = moments('bose', res, None, K)
Psi = lambda t: sum(Q(res[a]) / (t + Q(a) ** 2) for a in S)
D = lambda t: mp.fprod(t + Q(a) ** 2 for a in S)
worst = mp.mpf(0)
for e in range(K - 1):
    lhs = mp.quad(lambda y: wB(y) * y ** (2 * e) * Psi(y * y), BR)
    assert mom[e][2] == 0
    worst = max(worst, abs(lhs - (Q(mom[e][0]) + Q(mom[e][1]) * G)) / abs(lhs))
print(f"B1. K={K}: pi^2-coefficients of mu(t^e Psi), e=0..K-1:", [str(m[2]) for m in mom])
print("    max relative |quadrature - (c0 + cG G)| for e <= K-2:", mp.nstr(worst, 3))
worst = mp.mpf(0)
for t in (mp.mpf('0.3'), mp.mpf(2), mp.mpf(11), mp.mpf(70)):
    rhs = mp.quad(lambda s: s ** mp.mpf(-0.5) * (1 / D(t) - 1 / D(s)) / (s - t), [0, t / 2, t, 2 * t, 10 * t + 10, 1000, mp.inf]) / pi
    worst = max(worst, abs(rhs - Psi(t)) / Psi(t))
print("B2. max relative error of the integral representation of Psi at 4 points:", mp.nstr(worst, 3),
      "; min Psi on a grid of (0,10^4):", mp.nstr(min(Psi(mp.mpf(10) ** (k / 10.0)) for k in range(-30, 41)), 5))
c = sum(Q(r) for r in res.values())
print("B3. sum r_a a =", sum((res[a] * a for a in S), fmpq(0)), "; c = sum r_a =", mp.nstr(c, 8),
      "; int_0^inf (y^2 Psi(y^2) - c) dy =",
      mp.nstr(mp.quad(lambda y: -sum(Q(res[a]) * Q(a) ** 2 / (y * y + Q(a) ** 2) for a in S), [0, 1, 4, 16, 100, mp.inf]), 3),
      " [integrand evaluated in the cancellation-free form -sum r_a a^2/(y^2+a^2); pointwise check at y=3:",
      mp.nstr(9 * Psi(mp.mpf(9)) - c + sum(Q(res[a]) * Q(a) ** 2 / (9 + Q(a) ** 2) for a in S), 3), "]")
print("    signs of r_a:", ''.join('+' if r > 0 else '-' for r in res.values()), "; G-coefficient of mu(Psi) =", mom[0][1], "> 0")
