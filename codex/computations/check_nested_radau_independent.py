#!/usr/bin/env python3
"""Independent rational reconstruction of the nested Gaussian/Radau bounds.

Higher derivatives use Taylor-series division rather than the author's
rational-function differentiation recurrence. No numerical zeta values.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, prod


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(p, q, sign=1):
    out = [F(0)]*max(len(p), len(q))
    for j, a in enumerate(p):
        out[j] += a
    for j, b in enumerate(q):
        out[j] += sign*b
    return trim(out)


def mul(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for j, a in enumerate(p):
        for k, b in enumerate(q):
            out[j+k] += a*b
    return trim(out)


def value(p, x):
    return sum(c*x**j for j, c in enumerate(p))


def solve(matrix, rhs):
    a = [list(row)+[b] for row, b in zip(matrix, rhs)]
    for k in range(len(a)):
        assert a[k][k] != 0
        pivot = a[k][k]
        a[k] = [x/pivot for x in a[k]]
        for j in range(len(a)):
            if j != k:
                factor = a[j][k]
                a[j] = [x-factor*y for x, y in zip(a[j], a[k])]
    return [row[-1] for row in a]


bernoulli = [F(1)]
for n in range(1, 35):
    bernoulli.append(-sum(F(comb(n+1, k))*bernoulli[k]
                          for k in range(n))/F(n+1))


@lru_cache(None)
def moment(m, k):
    half = (F(2, 2**(2*k))-1)*bernoulli[2*k]
    return F((-1)**k, m*prod(range(1, 2*m, 2))) * half * prod(
        2*k+2*j+1 for j in range(m))


@lru_cache(None)
def gauss(m, n, shift):
    mu = [moment(m, j+shift) for j in range(2*n+1)]
    p = solve([[mu[j+k] for k in range(n)] for j in range(n)],
              [-mu[n+j] for j in range(n)])+[F(1)]
    assert all(sum(c*mu[j+k] for k, c in enumerate(p)) == 0
               for j in range(n))
    Q = [(-1)**(n+k)*c for k, c in enumerate(p)]
    assert all(c > 0 for c in Q)
    # Polynomial part of Q(z) times the Laurent series of the resolvent.
    P = [sum(Q[j]*(-1)**(j-k-1)*mu[j-k-1]
             for j in range(k+1, n+1)) for k in range(n)] or [F(0)]
    h = sum(c*mu[k] for k, c in enumerate(mul(p, p)))
    assert h > 0
    return p, P, Q, h


def divided_derivative(P, Q, order, z):
    def jet(poly):
        return [sum(c*comb(k, j)*z**(k-j)
                    for k, c in enumerate(poly) if k >= j)
                for j in range(order+1)]
    num, den = jet(P), jet(Q)
    coefficients = []
    for j in range(order+1):
        coefficients.append((num[j]-sum(den[k]*coefficients[j-k]
                                           for k in range(1, j+1)))/den[0])
    return (-1)**order*coefficients[-1]


wronskians = pairs = 0
for m in range(1, 5):
    assert moment(m, 0) == F(1, m)
    for shift in (0, 1):
        for n in range(7):
            _, P, Q, h = gauss(m, n, shift)
            _, P1, Q1, _ = gauss(m, n+1, shift)
            assert add(mul(P1, Q), mul(P, Q1), -1) == [h]
            wronskians += 1
    previous = None
    for n in range(7):
        p, P, Q, h = gauss(m, n, 0)
        pt, Pt, Qt, ht = gauss(m, n, 1)
        plusP = add([F(1, m)*c for c in Qt], Pt, -1)
        plusQ = [F(0)]+Qt
        endpoint = F(1, m)-Pt[0]/Qt[0]
        assert endpoint > 0
        endpoint_integral = sum(c*moment(m, j) for j, c in enumerate(mul(pt, pt)))
        assert endpoint == endpoint_integral/Qt[0]**2
        z = F((2*n+1)**2, 4)
        harmonic = sum((F(1, j**(2*m+1)) for j in range(1, n+1)), F(0))
        lower = harmonic+divided_derivative(P, Q, m-1, z)/2
        upper = harmonic+divided_derivative(plusP, plusQ, m-1, z)/2
        assert lower < upper
        if previous:
            assert previous[0] < lower < upper < previous[1]
        previous = lower, upper
        radius_lower = F(comb(2*n+m-1, m-1), 2)*h/(z**m*value(Q, z)**2)
        radius_upper = F(comb(2*n+m, m-1), 2)*ht/(z**(m+1)*value(Qt, z)**2)
        crude_lower = F(comb(2*n+m-1, m-1), 2)*moment(m, 2*n)/z**(2*n+m)
        crude_upper = F(comb(2*n+m, m-1), 2)*moment(m, 2*n+1)/z**(2*n+m+1)
        assert radius_lower <= crude_lower and radius_upper <= crude_upper
        assert upper-lower <= radius_lower+radius_upper
        if n == 1:
            s, cutoff = 2*m+1, 50
            finite = sum((F(1, j**s) for j in range(1, cutoff+1)), F(0))
            zeta_lower = finite+F(1, (s-1)*(cutoff+1)**(s-1))
            zeta_upper = finite+F(1, (s-1)*cutoff**(s-1))
            assert lower < zeta_lower < zeta_upper < upper
            assert zeta_upper-lower < radius_lower
            assert upper-zeta_lower < radius_upper
        pairs += 1

print(f"PASS: {wronskians} independently constructed Wronskians; {pairs} rational "
      "nested pairs including n=0; positive endpoint masses, both radii, and "
      "four enclosures verified against rational zeta-series intervals")
