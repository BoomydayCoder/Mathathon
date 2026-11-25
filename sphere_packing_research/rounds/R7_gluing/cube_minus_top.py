# Necessary conditions for N = D(M) minus {M}, M = product of odd primes <= y  (closed forms, relative to mass M).
import math
from sympy import primerange
def H(x): return -x*math.log(x)-(1-x)*math.log(1-x) if 0<x<1 else 0.0
for y in (11, 47, 200, 1000, 10000):
    P = list(primerange(3, y+1))
    theta = math.prod(1-1/p for p in P); mass = 1-theta
    Lam = sum(math.log(p)/(p-1) for p in P)
    # vertex condition at T = M/p (maximal in N): sum_{q!=p} log q/(q-1) >= H(x)/(2x), x = theta/((p-1)(1-theta))
    worst_v = min((Lam - math.log(p)/(p-1)) - H(theta/((p-1)*mass))/(2*theta/((p-1)*mass)) for p in P)
    # per-prime cut: matched mass along q = 1/q - theta/(q-1); x_q = (1/q)/mass
    worst_c = min(math.log(q)*(1/q - theta/(q-1)) - 0.5*mass*H((1/q)/mass) for q in P if (1/q)/mass < 1) if all((1/q)/mass<1 for q in P) else float('-inf')
    # singleton partition: I/M = sum_q log q (1/q - theta/(q-1));  (C/2)H(nu)/M = 0.5*[mass*log(mass) + E_mu[log(1/mu(T)); T != top]]
    Emu = sum((1/p)*math.log(p) + (1-1/p)*math.log(p/(p-1)) for p in P)      # full-cube entropy of mu
    ent = 0.5*(mass*math.log(mass) + Emu - theta*math.log(1/theta))
    I = sum(math.log(q)*(1/q - theta/(q-1)) for q in P)
    print(f"y={y:6d} #primes={len(P):5d} theta={theta:.4f} gain(cube-top)={1/mass:.4f} M/phi(M)={1/theta:.4f} | min vertex slack={worst_v:+.3f}  min per-prime cut slack={worst_c:+.5f}  singleton slack={I-ent:+.4f}")
