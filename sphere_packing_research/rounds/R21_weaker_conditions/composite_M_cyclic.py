# R21: two prime factors M = l1*l2.  Along the whole segment of class-weights c allowed by the prime-order conditions
# (between the extreme rays R2 [c_1 = 0] and R1 [c_a = 0]) test the two further necessary conditions for the cyclic
# subgroup H = C_M of the non-split torus (U = Fix(H), a = dim U, n = dim V):
#   (alpha)  a = 0  or  a >= n/|H|                       (glue index <= |H|^a)
#   (beta)   |det(1-g | U^perp)| >= (n/a)^a                (glue index <= coinvariants of the complement)
# Units: multiplicity of lambda in V|T' is S - c_lambda;  n = (M-1) S.
import math
from fractions import Fraction as Fr
def scan(l1, l2, steps=2000):
    M = l1 * l2; S = 1.0
    feas = []
    c1_max = S / l2                                  # ray R1
    for i in range(steps + 1):
        c1 = c1_max * i / steps
        cb = (S / l1 - c1) / (l2 - 1)                # class of order 2*l2  (l2-1 characters)
        ca = (S / l2 - c1) / (l1 - 1)                # class of order 2*l1
        cm = (S - c1 - (l1 - 1) * ca - (l2 - 1) * cb) / ((l1 - 1) * (l2 - 1))
        if min(ca, cb, cm) < -1e-12: continue
        n = (M - 1) * S; a = S - c1
        alpha = a >= n / M - 1e-12
        lhs = (S - ca) * math.log(l1) + (S - cb) * math.log(l2)
        beta = lhs >= a * math.log(n / a) - 1e-9
        gain = 2 * M * (c1 ** 2 + (l1 - 1) * ca ** 2 + (l2 - 1) * cb ** 2 + (l1 - 1) * (l2 - 1) * cm ** 2 - 1 / M)
        if alpha and beta: feas.append((c1 * M, gain))
    return feas
for (l1, l2) in [(3, 5), (3, 7), (5, 11), (3, 29), (7, 11), (11, 13)]:
    f = scan(l1, l2)
    print('M=%d*%d: feasible points on the segment (c_theta0 in units of S/M, T\'-gain/(p/(p-1))): %s'
          % (l1, l2, [(round(x, 4), round(g, 6)) for x, g in f][:5]), ' count', len(f))
