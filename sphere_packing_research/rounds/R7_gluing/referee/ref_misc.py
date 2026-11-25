import math, sys, itertools
import numpy as np
from scipy.optimize import brentq
from sympy import primerange, totient, divisors, factorint
from ref_lp import w_apostol, margin_y
def H2(x): return 0.0 if x <= 0 or x >= 1 else -x*math.log(x) - (1-x)*math.log(1-x)
# (1) x_hi(q): root of H2(x) = 2(1-x) log q/(q-1);  x_lo(q): root of H2(x) = 2 x log q (squarefree)
print("q, x_hi, q*x_hi, x_lo, q^2*x_lo, check x_hi<2/q")
tab = {}
for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 59):
    L = math.log(q)/(q-1)
    f = lambda x: H2(x) - 2*(1-x)*L
    xhi = brentq(f, 1e-12, 1-1e-12)
    g = lambda x: H2(x) - 2*x*math.log(q)
    xlo = brentq(g, 1e-15, 0.5) if q > 2 else float('nan')
    tab[q] = xhi
    print(f"{q:3d} {xhi:.4f} {q*xhi:.3f} {xlo:.5f} {q*q*xlo:.3f} {xhi < 2/q}")
print("sum x_hi {3,5,7} =", round(tab[3]+tab[5]+tab[7], 4), " {3,5,7,11} =", round(tab[3]+tab[5]+tab[7]+tab[11], 4),
      " {3,5,7,13} =", round(tab[3]+tab[5]+tab[7]+tab[13], 4), " {2,3,5} =", round(tab[2]+tab[3]+tab[5], 4), " {2,3,7} =", round(tab[2]+tab[3]+tab[7],4), " {2,5,7} =", round(tab[2]+tab[5]+tab[7],4))
# (2) brute-force per-prime cut slack for explicit sets (sanity of the vectorised filter in ref_even.py)
def perprime_slacks(N):
    M = math.lcm(*N); c = {n: int(totient(n)) for n in N}; C = sum(c.values()); out = {}
    for q, e in factorint(M).items():
        top = [n for n in N if factorint(n).get(q, 0) == e]; low = [n for n in N if n not in top]
        if not low: out[q] = None; continue
        Wc = sum(w_apostol(a, b) for a in top for b in low)
        out[q] = Wc - 0.5*C*H2(sum(c[n] for n in low)/C)
    return out
for N in ([1,3,5,7,21,35,105], [1,5,7,11,35,55,77,231,385,1155], [1,3,5,7,15,21,35], [4,12,20], [1,3,4,5,12,15,20]):
    print("per-prime cut slacks", N, {q: (None if s is None else round(s, 4)) for q, s in perprime_slacks(N).items()}, "LP margin", round(margin_y(N), 4))
# (3) cube minus top, M = product of odd primes <= y: per-prime cut slack (in units of M), my own closed form
for y in (11, 47, 200, 1000, 1500, 2000, 2500, 3000, 10000):
    P = list(primerange(3, y+1)); theta = math.prod(1-1/p for p in P); mass = 1-theta
    sl = []
    for q in P:
        Wq = math.log(q)*(1/q - theta/(q-1))          # edges T -> T+q with q not in T and T+q != top
        x = (1/q)/mass                                 # nu-fraction lacking q
        sl.append(Wq - 0.5*mass*H2(x))
    i = int(np.argmin(sl))
    print(f"y={y}: theta={theta:.4f} gain cube-minus-top={1/mass:.4f} M/phi(M)={1/theta:.4f} min per-prime cut slack={min(sl):+.5f} at q={P[i]}")
# exact brute-force check of that closed form for M=1155
M = 1155; N = [d for d in divisors(M) if d != M]
print("brute force M=1155 cube minus top per-prime slacks /M:", {q: round(s/M, 5) for q, s in perprime_slacks(N).items()})
P = [3,5,7,11]; theta = math.prod(1-1/p for p in P); mass = 1-theta
print("closed form:", {q: round(math.log(q)*(1/q - theta/(q-1)) - 0.5*mass*H2((1/q)/mass), 5) for q in P})
