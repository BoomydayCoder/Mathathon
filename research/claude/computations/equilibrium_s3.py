"""Verification of the weighted equilibrium problem with min E = (9/2) log 2.

Problem: minimise over probability measures rho on [0, inf)
    E(rho) = int V drho - iint log|x^2 - y^2| drho drho,
    V(x) = 2 pi x + int_0^1 log(x^2+u^2) du = 2 pi x + log(1+x^2) - 2 + 2x arctan(1/x).
Symmetrise: sigma(x) = rho(|x|) on [-b, b], mass 2; E = 1/2 [int V dsigma - iint log|x-y| dsigma dsigma].

Closed forms checked here (derivation summarised in comments):
  g(x) = sgn(x)(3pi/2 - arctan|x|) = int_0^inf w(a) x/(x^2+a^2) da,  w = 3 on (0,1), 2 on (1,inf).
  Each kernel x/(x^2+a^2) has the bounded-at-endpoints solution with Cauchy transform
      G_a(z) = z/(z^2+a^2) - R(z) a / (r_a (z^2+a^2)),  R(z)=sqrt(z^2-b^2) ~ z,  r_a = sqrt(a^2+b^2),
  density s a/(pi r_a (x^2+a^2)), s = sqrt(b^2-x^2), mass 1 - a/r_a.
  => sigma(x) = (1/pi) [3 artanh(s/b) - artanh(s/c)],  c = sqrt(1+b^2);  mass 3b + 1 - c = 2  =>  b = 3/4, c = 5/4.
  Off support (x > b):  2G(x) - V'(x) = 2[3 arcsin(b/x) - arctan(c/sqrt(x^2-b^2)) - pi] < 0.
  ell = 6b - 2c + 2 log(2/b) + 2 artanh(1/c) = 2 + 6 log 2.
  Virial (dilation) identity: int x V'(x) dsigma = mass^2 = 4  =>  int V dsigma = int log(1+x^2) dsigma.
  i G(it) = log(1+t) - 3 log t + 3 log(b+sqrt(t^2+b^2)) - log(c+sqrt(t^2+b^2)),
  int log(1+x^2) dsigma = 2U(0) + 2 int_0^1 iG(it) dt = (-2-ell) + 12 log 2 = 6 log 2 - 4.
  min E = (1/4) int V dsigma + ell/2 = (9/2) log 2.
"""
import mpmath as mp

mp.mp.dps = 40
TOL = mp.mpf(10) ** -25
b = mp.mpf(3) / 4
c = mp.sqrt(1 + b * b)
results = []


def check(name, val, ref, tol=TOL):
    err = abs(val - ref)
    ok = err < tol
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: value={mp.nstr(val, 30)}  ref={mp.nstr(ref, 30)}  err={mp.nstr(err, 3)}")


def check_true(name, cond, info=""):
    results.append(bool(cond))
    print(f"{'PASS' if cond else 'FAIL'}  {name} {info}")


def sig(x):
    x = abs(mp.mpf(x))
    if x >= b:
        return mp.mpf(0)
    s = mp.sqrt(b * b - x * x)
    return (3 * mp.log((b + s) / x) - mp.atanh(s / c)) / mp.pi  # artanh(s/b) = log((b+s)/x)


def sig_rep(x):  # integral representation over the charge parameter a
    s = mp.sqrt(b * b - x * x)
    h = lambda a: a / (mp.sqrt(a * a + b * b) * (x * x + a * a))
    return s / mp.pi * (3 * mp.quad(h, [0, x, 1]) + 2 * mp.quad(h, [1, mp.inf]))


def V(x):
    x = abs(mp.mpf(x))
    if x == 0:
        return mp.mpf(-2)
    return 2 * mp.pi * x + mp.log(1 + x * x) - 2 + 2 * x * mp.atan(1 / x)


def Vp(x):  # x > 0
    return 2 * mp.pi + 2 * mp.atan(1 / x)


def g(x):
    return mp.sign(x) * (3 * mp.pi / 2 - mp.atan(abs(x)))


def U(x):
    x = mp.mpf(x)
    pts = sorted(set([-b, -b / 2, 0, b / 2, b] + ([x] if -b < x < b else [])))
    return mp.quad(lambda y: (mp.log(abs(x - y)) if y != x else 0) * sig(y), pts)


print("=== (0) parameters ===")
check("V closed form = 2 pi x + int_0^1 log(x^2+u^2)du at x=0.4", V(0.4),
      2 * mp.pi * 0.4 + mp.quad(lambda u: mp.log(mp.mpf(0.4) ** 2 + u * u), [0, 1]))
check("mass condition 3b + 1 - sqrt(1+b^2) = 2 at b = 3/4", 3 * b + 1 - c, 2)
check("c = sqrt(1+b^2) = 5/4", c, mp.mpf(5) / 4)

print("=== (1) density ===")
for x in [mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf('0.6'), mp.mpf('0.74')]:
    check(f"sigma closed form = a-integral representation at x={x}", sig(x), sig_rep(x))
grid = [b * k / 400 for k in range(1, 400)]
check_true("sigma > 0 on (0,b) (grid of 399 pts; analytic: s/b > s/c and 3>1)", min(sig(x) for x in grid) > 0,
           f"min={mp.nstr(min(sig(x) for x in grid), 5)}")
mass = 2 * mp.quad(sig, [0, b / 2, b])
check("total mass of sigma = 2 (rho mass 1)", mass, 2)
for x in [mp.mpf('0.1'), mp.mpf('0.4'), mp.mpf('0.7')]:
    # PV int sigma(y)/(x-y) dy = int (sigma(y)-sigma(x))/(x-y) dy + sigma(x) log((b+x)/(b-x))
    pv = mp.quad(lambda y: (sig(y) - sig(x)) / (x - y), [-b, 0, x, b]) + sig(x) * mp.log((b + x) / (b - x))
    check(f"Hilbert transform PV int sigma/(x-y) = g(x) at x={x}", pv, g(x), mp.mpf(10) ** -20)

ell = V(b) - 2 * U(b)
for x in [mp.mpf(0), mp.mpf('0.1'), mp.mpf('0.3'), mp.mpf('0.5'), mp.mpf('0.7')]:
    check(f"E-L equality V - 2U = ell at x={x}", V(x) - 2 * U(x), ell, mp.mpf(10) ** -20)
for x in [mp.mpf('0.8'), mp.mpf(1), mp.mpf(2), mp.mpf(5), mp.mpf(20)]:
    gap = 2 * U(x) - V(x) + ell
    d_num = 2 * mp.quad(lambda y: sig(y) / (x - y), [-b, 0, b]) - Vp(x)
    d_cf = 2 * (3 * mp.asin(b / x) - mp.atan(c / mp.sqrt(x * x - b * b)) - mp.pi)
    check(f"(2U-V)' closed form off support at x={x}", d_num, d_cf, mp.mpf(10) ** -20)
    check_true(f"strict inequality 2U - V < -ell at x={x}", gap < 0 and d_cf < 0, f"2U-V+ell={mp.nstr(gap, 8)}")

print("=== (2) constants and energy ===")
check("ell = 2 + 6 log 2", ell, 2 + 6 * mp.log(2))
check("ell = 6b - 2c + 2log(2/b) + 2artanh(1/c)", ell, 6 * b - 2 * c + 2 * mp.log(2 / b) + 2 * mp.atanh(1 / c))
IV = 2 * mp.quad(lambda x: V(x) * sig(x), [0, b / 2, b])
check("int V dsigma = 6 log 2 - 4", IV, 6 * mp.log(2) - 4)
check("virial: int x V'(x) dsigma = 4", 2 * mp.quad(lambda x: x * Vp(x) * sig(x), [0, b / 2, b]), 4)
check("int V dsigma = int log(1+x^2) dsigma", IV, 2 * mp.quad(lambda x: mp.log(1 + x * x) * sig(x), [0, b / 2, b]))
iG = lambda t: mp.log(1 + t) - 3 * mp.log(t) + 3 * mp.log(b + mp.sqrt(t * t + b * b)) - mp.log(c + mp.sqrt(t * t + b * b))
for t in [mp.mpf('0.2'), mp.mpf(1), mp.mpf(3)]:
    check(f"iG(it) closed form at t={t}", iG(t), 2 * mp.quad(lambda y: t * sig(y) / (t * t + y * y), [0, b / 2, b]))
check("2 int_0^1 iG(it) dt = 12 log 2", 2 * mp.quad(iG, [0, 1]), 12 * mp.log(2))
E = IV / 4 + ell / 2
check("min E = (1/4) int V dsigma + ell/2 = (9/2) log 2", E, mp.mpf(9) / 2 * mp.log(2))
# direct double-integral energy as an independent check of the formula E = 1/2[intV - iint log]
mp.mp.dps = 25  # cheaper nested quadrature for this cross-check
I2 = 2 * mp.quad(lambda x: sig(x) * U(x), [0, b / 2, b])
mp.mp.dps = 40
check("direct E = 1/2[int V dsigma - iint log|x-y|] (tol 1e-18)", (IV - I2) / 2, mp.mpf(9) / 2 * mp.log(2), mp.mpf(10) ** -18)
print("min E =", mp.nstr(E, 30), "  (9/2)log2 =", mp.nstr(mp.mpf(9) / 2 * mp.log(2), 30))
print("ALL PASS" if all(results) else f"FAILURES: {results.count(False)}")
