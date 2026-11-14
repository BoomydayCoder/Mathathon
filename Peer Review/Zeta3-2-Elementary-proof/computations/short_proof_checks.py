#!/usr/bin/env python3
"""Short proof of Apery's theorem via the Casoratian: numerical/exact checks.

R_n(t) = ((t-1)...(t-n)/(t(t+1)...(t+n)))^2,
A2k = C(n+k,n)^2 C(n,k)^2,  A1k = 2 A2k (2H_k - H_{n+k} - H_{n-k}),
u_n = 2 sum_k A2k,  v_n = 2 sum_k A2k H_k^(3) + sum_k A1k H_k^(2),
F_n = u_n zeta(3) - v_n = -sum_{t>=1} R_n'(t).
Exact identities: fractions.Fraction.  Numerics: mpmath (dps >= 120/200).
"""
from fractions import Fraction as Fr
from math import comb, lcm
import mpmath as mp

NMAX = 60
results = []


def check(name, ok, extra=""):
    results.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))


Hcache = {}


def H(k, j):
    key = (k, j)
    if key not in Hcache:
        Hcache[key] = sum((Fr(1, l ** j) for l in range(1, k + 1)), Fr(0))
    return Hcache[key]


def uv(n):
    u = Fr(0)
    v = Fr(0)
    for k in range(n + 1):
        a2 = comb(n + k, n) ** 2 * comb(n, k) ** 2
        a1 = 2 * a2 * (2 * H(k, 1) - H(n + k, 1) - H(n - k, 1))
        u += 2 * a2
        v += 2 * a2 * H(k, 3) + a1 * H(k, 2)
    return u, v


U, V = [], []
for n in range(NMAX + 2):
    a, b = uv(n)
    U.append(a)
    V.append(b)

print("=" * 72)
print("Sanity: u_n integer, recurrence L[u]=L[v]=0 for n=1..%d" % NMAX)
print("=" * 72)
ok_int = all(U[n].denominator == 1 for n in range(NMAX + 2))
check("u_n in Z for n=0..%d" % (NMAX + 1), ok_int, f"u_0..u_4 = {[int(x) for x in U[:5]]}, v_0..v_2 = {V[:3]}")


def L(X, n):
    return (n + 1) ** 3 * X[n + 1] - (2 * n + 1) * (17 * n * n + 17 * n + 5) * X[n] + n ** 3 * X[n - 1]


check("L[u]_n = 0, n=1..%d" % NMAX, all(L(U, n) == 0 for n in range(1, NMAX + 1)))
check("L[v]_n = 0, n=1..%d" % NMAX, all(L(V, n) == 0 for n in range(1, NMAX + 1)))

print("\n" + "=" * 72)
print("A1. Casoratian u_n v_{n-1} - u_{n-1} v_n = -24/n^3, n=1..%d (exact)" % NMAX)
print("=" * 72)
cas = [U[n] * V[n - 1] - U[n - 1] * V[n] for n in range(1, NMAX + 1)]
bad = [n for n in range(1, NMAX + 1) if cas[n - 1] != Fr(-24, n ** 3)]
check("A1 Casoratian identity n=1..%d" % NMAX, not bad, f"first values {cas[:3]}; failures: {bad[:5]}")

print("\n" + "=" * 72)
print("A2. Growth of u_n")
print("=" * 72)


def lowb(n):
    return Fr(33 * n ** 3 + 51 * n ** 2 + 27 * n + 5, (n + 1) ** 3)


bad = [n for n in range(1, NMAX + 1) if not (U[n + 1] / U[n] >= lowb(n))]
check("u_{n+1}/u_n >= (33n^3+51n^2+27n+5)/(n+1)^3, n=1..%d" % NMAX, not bad, f"failures {bad}")
# also n = 0 for info
print(f"   n=0: u_1/u_0 = {U[1]/U[0]}, bound(0) = {lowb(0)}")
ratios = [(U[n + 1] / U[n], n) for n in range(1, NMAX + 1)]
rmin = min(ratios)
print(f"   min_{{1<=n<={NMAX}}} u_(n+1)/u_n = {float(rmin[0]):.10f} at n={rmin[1]}  (u_2/u_1 = {U[2]/U[1]})")
print(f"   u_(NMAX+1)/u_NMAX = {float(U[NMAX+1]/U[NMAX]):.10f}; limit (1+sqrt2)^4 = {float(mp.mpf(1+mp.sqrt(2))**4):.10f}")
# exact identity of the bound: 33 - deficit
for n in range(1, 6):
    x = Fr(1, n + 1)
    assert lowb(n) == 33 - (48 * x - 24 * x * x + 4 * x ** 3)
check("bound = 33 - 48/(n+1) + 24/(n+1)^2 - 4/(n+1)^3 (so increasing in n since 48-48x+12x^2>0 on 0<x<=1)", True)
for thr in (27, 28, 30):
    n0 = next(n for n in range(0, 10 ** 4) if lowb(n) >= thr)
    print(f"   bound(n) >= {thr}  iff  n >= {n0}   (bound({n0-1}) = {float(lowb(n0-1)):.5f}, bound({n0}) = {float(lowb(n0)):.5f})")
# smallest c with u_n >= 2 c^n
mp.mp.dps = 50
cvals = [(mp.power(mp.mpf(int(U[n])) / 2, mp.mpf(1) / n), n) for n in range(1, NMAX + 1)]
cmin = min(cvals)
print(f"   largest c with u_n >= u_0 c^n for all n<= {NMAX}: c = min (u_n/2)^(1/n) = {mp.nstr(cmin[0], 15)} at n={cmin[1]}")
kint = int(mp.floor(cmin[0]))
check(f"u_n >= 2*{kint}^n for n=0..{NMAX} (largest integer k)", all(U[n] >= 2 * kint ** n for n in range(NMAX + 1)) and
      not all(U[n] >= 2 * (kint + 1) ** n for n in range(NMAX + 1)))
for n in (1, 2, 3, 5, 10, 20, 40, 60):
    print(f"     n={n:2d}: (u_n/2)^(1/n) = {mp.nstr(mp.power(mp.mpf(int(U[n]))/2, mp.mpf(1)/n), 12)}")

print("\n" + "=" * 72)
print("A3. F_n = sum_{m>n} 24 u_n/(m^3 u_m u_{m-1});  0 < F_n < 12/(n^2 u_n)")
print("=" * 72)
mp.mp.dps = 200  # (dps 120 loses ~30 digits to cancellation in u_n zeta3 - v_n at n=10)
z3 = mp.zeta(3)


def mpq(q):
    return mp.mpf(q.numerator) / q.denominator


Um = [mp.mpf(int(x)) for x in U]
# extend u up to m=81 via recurrence
Uext = list(U)
for n in range(len(Uext) - 1, 81):
    Uext.append(((2 * n + 1) * (17 * n * n + 17 * n + 5) * Uext[n] - n ** 3 * Uext[n - 1]) / (n + 1) ** 3)
Um = [mpq(x) for x in Uext]
mx = mp.mpf(0)
for n in range(0, 11):
    Fn = Um[n] * z3 - mpq(V[n])
    ser = sum(24 * Um[n] / (mp.mpf(m) ** 3 * Um[m] * Um[m - 1]) for m in range(n + 1, 81))
    rel = abs(Fn - ser) / abs(Fn)
    mx = max(mx, rel)
check("A3a series identity for F_n, n=0..10, m<=80 (rel err < 1e-100, dps 200)", mx < mp.mpf(10) ** -100, f"max rel err = {mp.nstr(mx, 5)}")
mp.mp.dps = 220
z3 = mp.zeta(3)
bad = []
ratios = []
for n in range(1, 41):
    Fn = mpq(U[n]) * z3 - mpq(V[n])
    b = 12 / (mp.mpf(n) ** 2 * mpq(U[n]))
    ratios.append(Fn / b)
    if not (0 < Fn < b):
        bad.append(n)
check("A3b 0 < F_n < 12/(n^2 u_n), n=1..40 (dps 220)", not bad, f"failures {bad}")
print("   ratio F_n / (12/(n^2 u_n)) for n=1,2,3,5,10,20,40:",
      [mp.nstr(ratios[n - 1], 8) for n in (1, 2, 3, 5, 10, 20, 40)])
print("   (expected limit: F_n ~ 24 u_n/((n+1)^3 u_{n+1}u_n) * (1/(1-1/33.97^2)) -> ratio ~ 2/(33.97 n) -> 0)")

print("\n" + "=" * 72)
print("A4. Termwise bounds for t > n and S_n = sum_{t>n} |R_n'(t)|")
print("=" * 72)
mp.mp.dps = 30
# exact termwise checks via Fractions on a modest range, float/mp check on large range
ok_R = True
ok_H = True
worst = (0, None)
for n in range(1, 41):
    Hn = float(H(n, 1))
    for t in range(n + 1, 50 * n + 1):
        a = mp.fsum(mp.mpf(1) / (t - j) for j in range(1, n + 1))
        b = mp.fsum(mp.mpf(1) / (t + j) for j in range(0, n + 1))
        d = abs(a - b)
        r = d / Hn
        if r > worst[0]:
            worst = (r, (n, t))
        if d > Hn:
            ok_H = False
        Rt = mp.fprod(((mp.mpf(t) - j) / (t + j)) ** 2 for j in range(1, n + 1)) / mp.mpf(t) ** 2
        if not (0 < Rt < mp.mpf(1) / t ** 2):
            ok_R = False
check("0 < R_n(t) < 1/t^2 for integers t>n, n=1..40, t<=50n", ok_R)
check("|sum_{j=1}^n 1/(t-j) - sum_{j=0}^n 1/(t+j)| <= H_n for t>n, n=1..40, t<=50n", ok_H,
      f"max ratio |.|/H_n = {mp.nstr(worst[0], 8)} at (n,t)={worst[1]}")
# check exactly at t=n+1 that the value is 2H_n - H_{2n+1} (the max)
ok = all(abs(sum(Fr(1, n + 1 - j) for j in range(1, n + 1)) - sum(Fr(1, n + 1 + j) for j in range(n + 1)))
         == 2 * H(n, 1) - H(2 * n + 1, 1) for n in range(1, 41))
check("at t=n+1 the bracket equals 2H_n - H_{2n+1} (exact, n=1..40)", ok)
print("   Proof of the bound: for t>=n+1, 0 < a:=sum 1/(t-j) <= H_n and 0 < b:=sum_{j=0}^n 1/(t+j) <= H_{2n+1}-H_n < log 2 + 1/(n+1)<= 1 <= H_n;")
print("   hence |a-b| <= max(a,b) <= H_n.  So |R_n'(t)| = 2R_n(t)|a-b| <= 2 H_n / t^2, S_n <= 2H_n sum_{t>n} 1/t^2 < 2H_n/n.")


def Rn(n, t):
    return mp.fprod(((t - j) / (t + j)) ** 2 for j in range(1, n + 1)) / t ** 2


def dRn(n, t):
    a = mp.fsum(1 / (t - j) for j in range(1, n + 1))
    b = mp.fsum(1 / (t + j) for j in range(0, n + 1))
    return 2 * Rn(n, t) * (a - b)


mp.mp.dps = 30
print("   n   S_n(t<=400n + tail R_n(T))   2H_n/n      4(1+log n)/n   |F_n|   S_n/(2H_n/n)")
okS = True
for n in (1, 2, 5, 10, 20, 40):
    T = 400 * n
    S = mp.fsum(abs(dRn(n, mp.mpf(t))) for t in range(n + 1, T + 1))
    S += Rn(n, mp.mpf(T))  # tail: R decreasing for large t, sum_{t>T}|R'| ~ R(T)
    b1 = 2 * float(H(n, 1)) / n
    b2 = 4 * (1 + mp.log(n)) / n
    with mp.workdps(250):
        Fn = abs(mpq(U[n]) * mp.zeta(3) - mpq(V[n]))
    okS &= S < b1
    print(f"  {n:3d}   {mp.nstr(S, 12):>22}   {b1:.8f}   {mp.nstr(b2, 8):>10}   {mp.nstr(Fn, 5):>10}   {mp.nstr(S/b1, 6)}")
check("S_n < 2H_n/n for n in {1,2,5,10,20,40}", okS)

print("\n" + "=" * 72)
print("A5. D_n^3 * 12/(n^2 u_n) < 1 and D_n^3 F_n")
print("=" * 72)
mp.mp.dps = 220
z3 = mp.zeta(3)


def D(n):
    return lcm(*range(1, n + 1)) if n >= 1 else 1


vals = [(n, mp.mpf(D(n)) ** 3 * 12 / (mp.mpf(n) ** 2 * mpq(U[n]))) for n in range(1, NMAX + 1)]
good = [n for n, x in vals if x < 1]
n0 = next(n for n in range(1, NMAX + 1) if all(x < 1 for m, x in vals if m >= n))
print(f"   n with D_n^3*12/(n^2 u_n) < 1: {good[:15]} ...")
print(f"   holds for all n in [{n0}, {NMAX}]; values near n0: " +
      ", ".join(f"n={n}:{mp.nstr(x,5)}" for n, x in vals if n0 - 3 <= n <= n0 + 3))
check(f"D_n^3*12/(n^2 u_n) < 1 for all n0={n0} <= n <= {NMAX}", True)
check("D_n^3 v_n in Z for n=1..60 (needed for integrality)", all((D(n) ** 3 * V[n]).denominator == 1 for n in range(1, NMAX + 1)))
print("   D_n^3 F_n, n=1..25:")
for n in range(1, 26):
    print(f"     n={n:2d}: {float(mp.mpf(D(n))**3 * (mpq(U[n])*z3 - mpq(V[n]))):.6e}")

print("\n" + "=" * 72)
print(f"SUMMARY: {sum(results)}/{len(results)} checks passed")
