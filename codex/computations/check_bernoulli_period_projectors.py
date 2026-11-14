"""Exact operator checks and independent integral checks for the period tower."""
from math import factorial, comb
from fractions import Fraction as Q
import sympy as S
import mpmath as mp

t = S.symbols('t')
Z = S.Function('Z')


def bplus(j):
    return S.Rational(1, 2) if j == 1 else S.bernoulli(j)


def theta(f):
    return S.expand(t*S.diff(f, t))


for m in range(1, 9):
    exact = m*S.polylog(m+1, S.exp(t))-t*S.polylog(m, S.exp(t))
    assert S.simplify(S.diff(exact, t, m)-t/(1-S.exp(-t))) == 0
    y = sum((m-k)*Z(m+1-k)*t**k/S.factorial(k) for k in range(m))
    y += sum(bplus(j)*t**(m+j)/S.factorial(m+j) for j in range(13))
    reflection = 2*sum((m-k)*Z(m+1-k)*t**k/S.factorial(k)
                       for k in range(m) if (k-m) % 2)
    reflection += t**(m+1)/S.factorial(m+1)
    assert S.expand(y-(-1)**m*y.subs(t, -t)-reflection) == 0
    projected = y
    for k in range(1, m):
        projected = S.expand(projected-theta(projected)/k)
    taylor_projection = sum((-t)**k/S.factorial(k)*S.diff(y, t, k)
                            for k in range(m))
    assert S.expand(projected-taylor_projection) == 0
    expected = m*Z(m+1)+(-1)**(m-1)/S.factorial(m-1)*sum(
        bplus(j)*t**(m+j)/(S.factorial(j)*(m+j)) for j in range(13))
    assert S.expand(projected-expected) == 0
    for N in range(31):
        multiplier = Q(1)
        for k in range(1, m):
            multiplier *= 1-Q(N, k)
        assert multiplier.denominator == 1
        target = 1 if N == 0 else (0 if N < m else (-1)**(m-1)*comb(N-1, m-1))
        assert multiplier == target

for r in range(1, 6):
    for a in range(20):
        multiplier = Q(1)
        for k in range(1, r):
            multiplier *= 1-Q(a, k)
        assert multiplier.denominator == 1
print('PASS exact derivative tower, reflection, Euler/Taylor projectors and integer multipliers')

mp.mp.dps = 45
integrals = 0
for m in (1, 2, 4, 6):
    for a in (mp.mpf('0.25'), mp.mpf(1), mp.mpf(3)):
        q = mp.exp(-a)
        # Product-variable form of the original (m+1)-fold integral.
        def integrand(u):
            if u == 0:
                return mp.mpf(0)  # An endpoint value does not affect the integral.
            z = q*u
            ratio = 1 if z == 1 else -mp.log(z)/(1-z)
            return (-mp.log(u))**(m-1)*ratio/factorial(m-1)
        base = mp.quad(integrand, [0, 1])
        polylog_base = (m*mp.polylog(m+1, q)+a*mp.polylog(m, q))/q
        assert mp.almosteq(base, polylog_base, rel_eps=mp.mpf('1e-37'))
        tail = mp.quad(lambda x: x**m/mp.expm1(x), [a, mp.inf])/factorial(m-1)
        polylog_tail = m*sum(a**k/factorial(k)*mp.polylog(m+1-k, q)
                            for k in range(m+1))
        assert tail > 0
        assert mp.almosteq(tail, polylog_tail, rel_eps=mp.mpf('1e-37'))
        initial = m*mp.zeta(m+1)-mp.quad(
            lambda x: x**m/mp.expm1(x) if x else (1 if m == 1 else 0), [0, a]
        )/factorial(m-1)
        assert mp.almosteq(tail, initial, rel_eps=mp.mpf('1e-37'))
        integrals += 1
print(f'PASS {integrals} independent base-period and positive-tail integral comparisons')
