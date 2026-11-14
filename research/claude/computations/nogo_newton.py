#!/usr/bin/env python3
"""nogo_newton.py -- exact p-adic Newton polygons of the monic polynomial F = Delta_K/lead (W=1, nodes 1..K, odd s)
versus Theorem A:  under (H1) p>s, p^2>2K; (H2) (s-1)(p-1)>2K-2; (H3) m=1 or s-2m-1>m(2m-3), m=floor(K/p);
(H4) p does not divide num H_d^(s), d<=m:   #{roots with v_p <= -s} >= n_p - m  (n_p=K-p+1; equality n_p, all =-s, if m=1).
Also prints the observed full pattern and cost_p = sum max(0,-v).  Exact computation (python-flint).
Usage: python3 nogo_newton.py s K [s K ...]"""
import sys, math
sys.path.insert(0, '/home/user/Mathathon/research/claude/computations')
from zeta7_lib import *
from fractions import Fraction
def newton_slopes(cs, p):
    pts = [(k, vp(c, p)) for k, c in enumerate(cs) if c != 0]
    hull = []
    for pt in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (y2 - y1) * (pt[0] - x1) >= (pt[1] - y1) * (x2 - x1): hull.pop()
            else: break
        hull.append(pt)
    return [(Fraction(-(y2 - y1), (x2 - x1)), x2 - x1) for (x1, y1), (x2, y2) in zip(hull, hull[1:])]
def Hnum(d, s):
    return sum((fmpq(1, v ** s) for v in range(1, d + 1)), fmpq(0)).p
def run(s, K):
    r = analyse(s, list(range(1, K + 1)), W=None, h=K)
    cs = r['cs'][:r['deg'] + 1]; lc = cs[-1]; cs = [c / lc for c in cs]
    assert abs(lc) == fmpq(math.factorial(K)) ** (s - 1), "lead != K!^(s-1)"
    gam = -log_abs_fmpq(content(cs)) / K ** 2; dl = -(r['logdet'] - log_abs_fmpq(lc)) / K ** 2
    print(f"\n=== s={s} K={K}: lead=+-K!^(s-1) OK; gamma=log lc(P_K)/K^2={gam:.5f}; delta=-log|F(zeta)|/K^2={dl:.5f}; total={gam-dl:+.5f}")
    print("  p   m  n_p  #(v<=-s)  cost  s*(n_p-m)  (s+2)n_p-2m  hyp  check   root valuations^multiplicity")
    proved = 0.0; viol = 0
    for p in primes_upto(K):
        if p * p <= 2 * K or p <= s: continue
        m = K // p; n_p = K - p + 1
        sl = newton_slopes(cs, p)
        cnt = sum(mu for v, mu in sl if v <= -s); cost = sum(-v * mu for v, mu in sl if v < 0)
        H2 = (s - 1) * (p - 1) > 2 * K - 2; H3 = (m == 1) or (s - 2 * m - 1 > m * (2 * m - 3))
        H4 = all(Hnum(d, s) % p != 0 for d in range(1, m + 1))
        hyp = H2 and H3 and H4
        ok = ''
        if hyp:
            good = cnt >= n_p - m and (m > 1 or (cnt == n_p and all(v >= -s for v, mu in sl)))
            ok = 'OK' if good else 'VIOLATION'; viol += (not good); proved += s * (n_p - m) * math.log(p)
        tag = ('A' if hyp else '-') + ('' if H4 else ' exc')
        print(f"{p:4d} {m:3d} {n_p:4d} {cnt:6d} {str(cost):>8s} {s*(n_p-m):8d} {(s+2)*n_p-2*m:9d}   {tag:5s} {ok:9s} " + ' '.join(f"{v}^{mu}" for v, mu in sl))
    print(f"  Theorem A proves log lc(P_K) >= {proved:.1f} = {proved/K**2:.4f} K^2   (true {gam:.4f} K^2); violations: {viol}")
a = sys.argv[1:]
for i in range(0, len(a), 2): run(int(a[i]), int(a[i + 1]))
