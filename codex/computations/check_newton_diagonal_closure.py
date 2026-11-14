"""Exact checks of the new Newton-kernel adjacent-diagonal closure.

The general proof is in newton-symmetric-square-review.md; finite checks
supplement the symbolic rational identity and check its boundary conventions.
"""
from math import comb
from fractions import Fraction
import sympy as S

n, k, c = S.symbols('n k c')
r = (n-k)/(n+k)
t_k_over_f2 = k**4/(n**2*(n+k)**2)
t_next_over_f2 = (n-k)**2/n**2
assert S.cancel(1+r*r-6*r-4*(t_k_over_f2-t_next_over_f2)) == 0
assert S.solve(2+c+2*(2-c), c) == [6]
middle = 6*(2*n+1)*(3*n*(n+1)+1)-(n+1)**3-n**3
assert S.expand(middle-(2*n+1)*(17*n*(n+1)+5)) == 0


def phi(n, k):
    return comb(n, k)*comb(n+k, k) if 0 <= k <= n else 0


def kernel(n, m):
    return sum(phi(n, k)*phi(m, k) for k in range(min(n, m)+1))


def tvalue(n, k):
    return Fraction(k**4*phi(n, k)**2, n*n*(n+k)**2)


for n in range(1, 41):
    assert tvalue(n, 0) == tvalue(n, n+1) == 0
    for k in range(n+1):
        f, g = phi(n, k), phi(n-1, k)
        assert f*f+g*g-6*f*g == 4*(tvalue(n, k)-tvalue(n, k+1))
    for k in range(1, n+2):
        binomial = comb(n-1, k-1)*comb(n+k-1, k-1) if k <= n else 0
        assert tvalue(n, k) == binomial**2
    assert kernel(n, n)+kernel(n-1, n-1) == 6*kernel(n, n-1)
    weights = [(comb(n, k)*comb(n+k-1, k))**2 for k in range(n+1)]
    assert 2*sum(k*k*w for k, w in enumerate(weights)) == n*n*sum(weights)
    assert kernel(n, n)+kernel(n-1, n-1) == 3*sum(weights)
    assert 2*kernel(n, n-1) == sum(weights)
    assert ((n+1)**3*kernel(n+1, n+1)
            == (2*n+1)*(17*n*(n+1)+5)*kernel(n, n)
               - n**3*kernel(n-1, n-1))

assert kernel(0, 0) == 1 and kernel(1, 1) == 5
for m in range(9):
    poly = [comb(m, j)**2 for j in range(m+1)]
    square = [sum(poly[j]*poly[d-j] for j in range(m+1)
                  if 0 <= d-j <= m) for d in range(2*m+1)]
    for n in range(16):
        coeff = sum(a*comb(n-d+2*m, 2*m) for d, a in enumerate(square) if d <= n)
        assert coeff == kernel(n, m)
print('PASS symbolic telescoper, forced coefficient6, and exact n=1..40 closure/recurrence')
print('PASS Legendre-square generating coefficients for m=0..8, n=0..15')
