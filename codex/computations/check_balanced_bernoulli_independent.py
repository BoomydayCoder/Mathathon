#!/usr/bin/env python3
"""Independent arbitrary-input checks of the balanced Bernoulli reduction.

Tests the parameterization on arbitrary V, not just solved balanced systems,
and tests the moment identity on arbitrary p with nonvanishing moments.
No constant-term conjecture is assumed.
"""
from math import factorial
import sympy as sp

x, z, t = sp.symbols("x z t")


def bernoulli(poly):
    return sum(c * sp.bernoulli(e[0], 0)
               for e, c in sp.Poly(sp.expand(poly), x).terms())


def mu(j, poly):
    return sum(
        c * (-1)**e[0] * sp.bernoulli(2 * e[0] + 2, 0)
        * sp.Rational(factorial(2 * e[0] + 2 * j + 1),
                      factorial(2 * j) * factorial(2 * e[0] + 2))
        for e, c in sp.Poly(sp.expand(poly), t).terms()
    )


def to_t(poly):
    terms = sp.Poly(sp.expand(poly), x).terms()
    assert all(e[0] % 2 == 0 for e, c in terms)
    return sum(c * (-t)**(e[0] // 2) for e, c in terms)


count = 0
nonzero_higher = 0
for m in range(1, 4):
    for n in range(1, 4):
        d, N, degree = m * n, 2 * m + 1, m * (n - 1)
        scale = sp.Rational((-1)**d, 2 * (d + m + 1))
        V = scale * z**degree + sum((-1)**k * (k + 2) * z**k
                                   for k in range(degree))
        Q = sp.expand((x * (x + 1))**N * V.subs(z, x * (x + 1)))
        R = sp.expand(Q - Q.subs(x, x - 1))
        assert sp.expand(Q - Q.subs(x, -x - 1)) == 0
        assert sp.expand(R - Q + Q.subs(x, -x)) == 0
        p = to_t(sp.cancel(R / x**N))
        assert sp.Poly(p, t).degree() == d and sp.Poly(p, t).LC() == 1
        assert p.subs(t, 0) == 2 * V.subs(z, 0)
        for j in range(1, m + 1):
            assert mu(j, t**(m - j) * p) == 0
            assert bernoulli(sp.diff(x*x*R, x, 2*j - 1)) == -bernoulli(
                sp.diff((2*x + 1)*Q, x, 2*j - 1))
            if mu(j, t**(m - j + 1) * p) != 0:
                nonzero_higher += 1
        # Independently check the coefficient identity away from its zero set.
        arbitrary = 3 + sum((-1)**k * (k + 1) * t**k for k in range(1, d + 1))
        arbitrary_R = sp.expand(x**N * arbitrary.subs(t, -x*x))
        for j in range(1, m + 1):
            for k in range(3):
                right = sp.Rational((-1)**(m - j + k), factorial(2*j))
                right *= bernoulli(sp.diff(x**(2*k)*arbitrary_R, x, 2*j - 1))
                assert mu(j, t**(m - j + k)*arbitrary) == right
                count += 1
assert nonzero_higher > 0
print(f"PASS: nine arbitrary V parameterizations; {count} arbitrary-p moment "
      f"identities; {nonzero_higher} nonzero higher constraints retained")
