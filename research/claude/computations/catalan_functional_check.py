#!/usr/bin/env python3
"""
catalan_functional_check.py -- numerical verification (50 digits working precision) of the
integral representations behind catalan_lib.py.  CHECKED NUMERICALLY ONLY (the proofs are the
classical Hermite / Abel-Plana formulas plus s-1 integrations by parts; see the report).

For kernel in {bose, sinh}, several s and alpha:
  (M)  int_0^inf y^(2e) w(y) dy          == mu(t^e)                 (rational, Bernoulli)
  (P)  int_0^inf w(y)/(y^2+a^2) dy       == const + coefX * X0      (a = alpha+j)
  (W)  w(y) > 0 on a grid; series form   w = c * y^s sum_l (c_l)^(s-1) e^{-c_l y}  with positive terms
"""
import sys
import mpmath as mp
from flint import fmpq
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from catalan_lib import mu_monomial, pole_values

mp.mp.dps = 50
TOL = mp.mpf(10) ** -32
fails = 0


def f(kernel, y):
    return 1 / mp.expm1(2 * mp.pi * y) if kernel == 'bose' else 1 / mp.sinh(mp.pi * y)


def w(kernel, s, y):
    y = mp.mpf(y)
    if s == 1:
        d = f(kernel, y)
    else:
        d = mp.diff(lambda u: f(kernel, u), y, s - 1)
    c = 2 if kernel == 'bose' else 1
    return c * (-1) ** (s - 1) * y ** s * d / mp.factorial(s - 1)


def w_series(kernel, s, y):
    """closed positive series: bose: 2 y^s sum_{l>=1} (2 pi l)^(s-1) e^{-2 pi l y}/(s-1)!;
       sinh: 2 y^s sum_{l odd} (pi l)^(s-1) e^{-pi l y}/(s-1)!"""
    y = mp.mpf(y)
    if kernel == 'bose':
        return 2 * y ** s / mp.factorial(s - 1) * mp.nsum(lambda l: (2 * mp.pi * l) ** (s - 1) * mp.exp(-2 * mp.pi * l * y), [1, mp.inf])
    return 2 * y ** s / mp.factorial(s - 1) * mp.nsum(lambda l: (mp.pi * (2 * l + 1)) ** (s - 1) * mp.exp(-mp.pi * (2 * l + 1) * y), [0, mp.inf])


def X0(kernel, s, alpha):
    alpha = mp.mpf(int(alpha.p)) / int(alpha.q)
    if kernel == 'bose':
        return mp.zeta(s, alpha)
    if s == 1:
        return (mp.digamma((alpha + 1) / 2) - mp.digamma(alpha / 2)) / 2
    return (mp.zeta(s, alpha / 2) - mp.zeta(s, (alpha + 1) / 2)) / mp.mpf(2) ** s


def q2m(q):
    return mp.mpf(int(q.p)) / mp.mpf(int(q.q))


def report(name, val, ref, tol=None):
    global fails
    err = abs(val - ref)
    ok = err < (tol or TOL) * max(1, abs(ref))
    fails += (not ok)
    print(f"{'PASS' if ok else 'FAIL'} {name}: {mp.nstr(val, 34)}  ref {mp.nstr(ref, 34)}  err {mp.nstr(err, 3)}")


cases = [('sinh', 2, fmpq(1, 2)), ('sinh', 1, fmpq(1)), ('sinh', 1, fmpq(1, 2)), ('sinh', 4, fmpq(1, 2)),
         ('sinh', 3, fmpq(1)), ('sinh', 2, fmpq(1, 3)), ('bose', 2, fmpq(1)), ('bose', 2, fmpq(1, 4)),
         ('bose', 2, fmpq(1, 3)), ('bose', 3, fmpq(1))]
print("constants: 4*Catalan =", mp.nstr(4 * mp.catalan, 30), " eta(2,1/2) =", mp.nstr(X0('sinh', 2, fmpq(1, 2)), 30))
report("eta(2,1/2) = 4 G", X0('sinh', 2, fmpq(1, 2)), 4 * mp.catalan)
report("eta(1,1) = log 2", X0('sinh', 1, fmpq(1)), mp.log(2))
report("eta(1,1/2) = pi/2", X0('sinh', 1, fmpq(1, 2)), mp.pi / 2)
report("zeta(2,1/4) = pi^2 + 8 G", X0('bose', 2, fmpq(1, 4)), mp.pi ** 2 + 8 * mp.catalan)
for (kernel, s, alpha) in cases:
    print(f"--- kernel={kernel} s={s} alpha={alpha} ---")
    x0 = X0(kernel, s, alpha)
    brk = [0, mp.mpf(1) / 4, 1, 3, 8, 20, 60]
    for e in (0, 1, 2, 5):
        val = mp.quad(lambda y: y ** (2 * e) * w(kernel, s, y), brk)
        report(f"(M) e={e}", val, q2m(mu_monomial(e, s, kernel)))
    pv = pole_values(s, kernel, alpha, 6)
    for j in (0, 1, 2, 5):
        a = q2m(alpha + j)
        val = mp.quad(lambda y: w(kernel, s, y) / (y * y + a * a), brk)
        report(f"(P) a={alpha + j}", val, q2m(pv[j][0]) + q2m(pv[j][1]) * x0)
    for y in (mp.mpf('0.01'), mp.mpf('0.3'), 1, 4):
        report(f"(W) series form at y={y}", w(kernel, s, y), w_series(kernel, s, y))
    grid_ok = all(w(kernel, s, mp.mpf(k) / 50) > 0 for k in range(1, 400))
    fails += (not grid_ok)
    print(f"{'PASS' if grid_ok else 'FAIL'} (W) w>0 on grid (0,8]")
# ---- functional II for Catalan ('cosh') and the chi_{-3} functional ('chi3'): explicit weights ----
w_cosh = lambda u: mp.pi * u / 4 * mp.sinh(mp.pi * u) / mp.cosh(mp.pi * u) ** 2
kap3 = lambda u: 1 / (3 + 4 * mp.sinh(mp.pi * u) ** 2)
w_chi3 = lambda u: 2 / (3 * mp.sqrt(3)) * u * 8 * mp.pi * mp.sinh(mp.pi * u) * mp.cosh(mp.pi * u) / (3 + 4 * mp.sinh(mp.pi * u) ** 2) ** 2   # = -(2/(3 sqrt3)) u kap3'(u)
L3 = (mp.zeta(2, mp.mpf(1) / 3) - mp.zeta(2, mp.mpf(2) / 3)) / 9
for (kernel, wf, x0) in (('cosh', w_cosh, mp.catalan), ('chi3', w_chi3, L3)):
    print(f"--- kernel={kernel} s=2 nodes n=1,2,... ---")
    brk = [0, mp.mpf(1) / 4, 1, 3, 8, 20, 60]
    for e in (0, 1, 2, 5):
        report(f"(M) e={e}", mp.quad(lambda y: y ** (2 * e) * wf(y), brk), q2m(mu_monomial(e, 2, kernel)))
    pv = pole_values(2, kernel, fmpq(1), 6)
    for j in (0, 1, 2, 5):
        report(f"(P) n={j + 1}", mp.quad(lambda y: wf(y) / (y * y + (j + 1) ** 2), brk), q2m(pv[j][0]) + q2m(pv[j][1]) * x0)
    grid_ok = all(wf(mp.mpf(k) / 50) > 0 for k in range(1, 400))
    fails += (not grid_ok)
    print(f"{'PASS' if grid_ok else 'FAIL'} (W) w>0 on grid (0,8]")
# ---- end-to-end: Hankel determinant from quadrature of the Gram matrix vs exact Delta_K(X0) from catalan_lib ----
from catalan_lib import analyse
w_sinh2 = lambda u: mp.pi * u * u * mp.cosh(mp.pi * u) / mp.sinh(mp.pi * u) ** 2
for (kernel, alpha, wf, x0) in (('cosh', fmpq(1), w_cosh, mp.catalan), ('sinh', fmpq(1, 2), w_sinh2, 4 * mp.catalan), ('chi3', fmpq(1), w_chi3, L3)):
    K = 4
    r = analyse(2, kernel, alpha, K)
    exact = sum(q2m(c) * x0 ** k for k, c in enumerate(r['cs']))
    al = q2m(alpha)
    D = lambda y: mp.fprod(y * y + (al + j) ** 2 for j in range(K))
    mom = [mp.quad(lambda y: y ** (2 * e) * wf(y) / D(y), [0, mp.mpf(1) / 4, 1, 3, 8, 20, 60]) for e in range(2 * K - 1)]
    num = mp.det(mp.matrix([[mom[i + j] for j in range(K)] for i in range(K)]))
    report(f"(D) K=4 Gram determinant by quadrature vs exact polynomial at X0, kernel={kernel} (tol 1e-25: the determinant amplifies quadrature error)", num / exact, mp.mpf(1), tol=mp.mpf(10) ** -25)
print("ALL PASS" if fails == 0 else f"{fails} FAILURES")
