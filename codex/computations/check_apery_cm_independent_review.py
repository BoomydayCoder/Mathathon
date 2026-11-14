#!/usr/bin/env python3
"""Independent exact checks supplementing the proof audit, not its substitute."""
from math import comb
from fractions import Fraction as F
import sympy as S


def summand(n,k):
    return comb(n,k)**2*comb(2*k,k) if 0 <= k <= n else 0


def certificate(n,k):
    return F(summand(n+1,k)*k**3*(3*k-4*n-4),(n+1)**2)


count = 0
for n in range(1,31):
    for k in range(n+2):
        lhs = ((n+1)**2*summand(n+1,k)
               -(10*n*n+10*n+3)*summand(n,k)
               +9*n*n*summand(n-1,k))
        assert lhs == certificate(n,k+1)-certificate(n,k)
        count += 1
assert count == 525

x = S.symbols('x')
M = (x-S.Rational(1,9))/(x-1)
h = 1/(1-x)
P = lambda z: (27*z*z-20*z+1)/(z*(z-1)*(9*z-1))
Q = lambda z: (9*z-3)/(z*(z-1)*(9*z-1))
Ps = P(M)*S.diff(M,x)-S.diff(M,x,2)/S.diff(M,x)
Qs = Q(M)*S.diff(M,x)**2
assert S.factor(M.subs(x,M)-x) == 0
assert S.factor(Ps-2*S.diff(h,x)/h-P(x)) == 0
assert S.factor(Qs-Ps*S.diff(h,x)/h+2*(S.diff(h,x)/h)**2
                -S.diff(h,x,2)/h-Q(x)) == 0

for p in [5,7]:
    bs = [sum(summand(n,k) for k in range(n+1)) % p for n in range(p)]
    B = S.Poly(sum(c*x**i for i,c in enumerate(bs)),x,modulus=p)
    D = S.Poly(9*x*x-18*x+1,x,modulus=p)
    assert S.gcd(B,D).degree() == 0
    assert S.legendre_symbol(-6,p) == 1
    assert B.eval(1) % p == 1
    print('Small prime',p,'critical remainder:',B.rem(D).as_expr())

assert S.factorint(4834944**2-4*14670139392) == {2:19,3:6,13:2,19:2}
print('PASS 525 direct finite identities including boundaries; rational symmetry; small-prime checks')
