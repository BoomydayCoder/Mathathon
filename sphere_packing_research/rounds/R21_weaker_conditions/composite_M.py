# R21: SL_2(p), p = 1 mod 4, M = (p+1)/2 squarefree composite.  V = sum_theta w_theta X_theta over odd characters theta of the
# non-split torus T' = C_{2M} (X_theta0 = psi + psi').  c_theta := coefficient of theta in  -chi_V|_{T' minus centre}
# (c constant on Galois classes = orders 2d, d | M).  Necessary condition (prime-order elements of T' regular; "free" is
# impossible when M is composite): push-forward of c to the characters of C_l is uniform, for every prime l | M.
# Gain carried by T':  G_T' = (2pM/(p-1)) (sum c^2 / S^2 - 1/M),  S = sum c.   We enumerate the extreme rays of the cone of
# admissible c (exact rational arithmetic) and maximise the convex function sum c^2 / S^2 over them.
import itertools, sys
from fractions import Fraction as Fr
from sympy import factorint, totient, divisors, Matrix, isprime

def analyse(p):
    M = (p + 1) // 2
    fac = factorint(M)
    assert all(e == 1 for e in fac.values()) and len(fac) >= 2
    primes = sorted(fac)
    ds = divisors(M)                                   # classes of odd characters: order 2d
    size = {d: int(totient(d)) for d in ds}             # number of characters in class d
    # uniformity at l:  sum_{d: l not | d} size_d c_d  =  (1/(l-1)) sum_{d: l | d} size_d c_d
    A = [[Fr(size[d]) if d % l else Fr(-size[d], l - 1) for d in ds] for l in primes]
    k = len(primes); best = None; rays = []
    for supp in itertools.combinations(range(len(ds)), k + 1):
        sub = Matrix([[A[i][j] for j in supp] for i in range(k)])
        ns = sub.nullspace()
        if len(ns) != 1: continue
        v = [Fr(int(x.p), int(x.q)) for x in ns[0]]
        if all(x <= 0 for x in v): v = [-x for x in v]
        if not all(x >= 0 for x in v) or not any(v): continue
        c = [Fr(0)] * len(ds)
        for j, x in zip(supp, v): c[j] = x
        S = sum(size[d] * c[i] for i, d in enumerate(ds))
        Q = sum(size[d] * c[i] ** 2 for i, d in enumerate(ds))
        g = Fr(2 * p * M, p - 1) * (Q / S ** 2 - Fr(1, M))
        rays.append((g, {d: c[i] / S for i, d in enumerate(ds) if c[i]}))
    rays.sort(key=lambda t: -t[0])
    return M, primes, rays

for p in [29, 41, 109, 173, 281, 389, 409, 1009, 2729]:
    M = (p + 1) // 2
    f = factorint(M)
    if not (isprime(p) and p % 4 == 1 and len(f) >= 2 and all(e == 1 for e in f.values())): continue
    M, primes, rays = analyse(p)
    g, c = rays[0]
    line = 'p=%d  M=%d=%s : max T\'-gain over %d extreme rays = %s = %.4f' % (p, M, '*'.join(map(str, primes)), len(rays), g, float(g))
    if len(primes) == 2:
        l1, l2 = primes
        formula = Fr(2 * p, p - 1) * Fr(l1 - 1, l2 - 1)
        line += '   [formula 2p(l1-1)/((p-1)(l2-1)) = %.4f, equal: %s]' % (float(formula), formula == g)
    print(line)
    print('      total gain <= 1 (centre) + (p+1)/(p-1) (unipotent) + this = %.3f   versus  p(p+1)/(p-1) = %.1f if M were prime' % (1 + (p + 1) / (p - 1) + float(g), p * (p + 1) / (p - 1)))
