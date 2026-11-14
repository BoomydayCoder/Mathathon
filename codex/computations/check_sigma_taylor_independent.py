#!/usr/bin/env python3
"""Independent rational audit of sigma jets, Casoratians and 20 certificates.

Uses the finite Newton quotient directly, rather than constructing the full
polynomial in t or generating its numerator from a recurrence. Standard library.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd, prod
from pathlib import Path
import json

R = 4


def constant(c):
    return [F(c)]+[F(0)]*R


def affine(c, d):
    return [F(c), F(d)]+[F(0)]*(R-1)


def add(a, b):
    return [x+y for x, y in zip(a, b)]


def scale(a, c):
    return [c*x for x in a]


def multiply(a, b):
    return [sum(a[j]*b[k-j] for j in range(k+1)) for k in range(R+1)]


def quotient(a, b):
    out = []
    for k in range(R+1):
        out.append((a[k]-sum(b[j]*out[k-j] for j in range(1, k+1)))/b[0])
    return out


def build(n):
    weights = [constant(1)]
    products = [constant(1)]
    for k in range(1, n+1):
        weights.append(quotient(weights[-1], affine(1, F(-1, k*k))))
        products.append(multiply(products[-1], affine(k*k, -1)))
    A, N = constant(0), constant(0)
    for k in range(n+1):
        term = scale(weights[k], comb(n, k)**2*comb(n+k, k)**2)
        A = add(A, term)
        inner = constant(0)
        for m in range(1, k+1):
            inner = add(inner, scale(products[m-1],
                F((-1)**(m-1)*factorial(n-m), m*factorial(n+m))))
        N = add(N, multiply(term, inner))
    tail = [sum((F(1, h**(2*r+3)) for h in range(1, n+1)), F(0))
            for r in range(R+1)]
    b = add(multiply(A, tail), scale(N, F(1, 2)))
    return A, b, quotient(b, A)


def alpha(n):
    k = n+1
    return scale(quotient(multiply(affine(k*k, -1), affine(k*k, -1)),
                          affine(3*k*k, -1)), F(k, 2))


data = [build(n) for n in range(14)]
assert data[0] == (constant(1), constant(0), constant(0))
assert data[1][2] == [F(1)+F(1, 5**(r+1)) for r in range(R+1)]
for n in range(1, 13):
    be = add(scale(affine(3*n*(n+1)+1, -1), 2*n+1),
             scale(add(alpha(n), alpha(n-1)), -1))
    for component in (0, 1):
        residual = add(multiply(alpha(n), data[n+1][component]),
                       multiply(alpha(n-1), data[n-1][component]))
        residual = add(residual, scale(multiply(be, data[n][component]), -1))
        assert residual == constant(0)
    W = add(multiply(data[n][1], data[n-1][0]),
            scale(multiply(data[n-1][1], data[n][0]), -1))
    assert W == quotient(constant(1), alpha(n-1))

# Independently reconstruct the exact Euler--Maclaurin enclosures in the JSON.
B = [F(1)]
for k in range(1, 121):
    B.append(-sum(comb(k+1, j)*B[j] for j in range(k))/F(k+1))


def zeta_interval(s, N=80, L=60):
    center = sum((F(1, j**s) for j in range(1, N)), F(0))
    center += F(1, (s-1)*N**(s-1))+F(1, 2*N**s)
    for k in range(1, L+1):
        center += B[2*k]*F(prod(range(s, s+2*k-1)),
                            factorial(2*k)*N**(s+2*k-1))
    sup_bound = sum(comb(2*L, j)*abs(B[j]) for j in range(2*L+1))
    radius = sup_bound*F(prod(range(s, s+2*L-1)),
                         factorial(2*L)*N**(s+2*L-1))
    return center-radius, center+radius


intervals = {r: zeta_interval(2*r+3) for r in (1, 2)}
rows = json.loads(Path(__file__).with_name('sigma_taylor_zeta_experiment.json')
                  .read_text())['rows']
assert len(rows) == 20
small = []
for row in rows:
    n, r = row['n'], row['r']
    A, b, q = data[n]
    q = q[r]
    assert q == F(row['q'])
    assert A[0] == F(row['A0']) and A[0].denominator == 1
    bad = factorial(n)*abs(A[1].numerator*b[0].numerator)
    # The actual finite experiment computes jets through order two.
    for x in A[:3]+b[:3]:
        bad *= x.denominator
    good = A[0].numerator
    while gcd(good, bad) > 1:
        good //= gcd(good, bad)
    assert good == int(row['good_factor'])
    forced = good**(r+1)
    assert forced == int(row['forced_component'])
    denominator = q.denominator
    assert denominator % forced == 0
    assert gcd(denominator//forced, good) == 1
    assert len(str(denominator)) == row['denominator_digits']
    lo, hi = intervals[r]
    if r == 1:
        assert q > hi and row['error_sign'] == 'above'
        elo, ehi = q-hi, q-lo
    else:
        assert q < lo and row['error_sign'] == 'below'
        elo, ehi = lo-q, hi-q
    assert elo == F(row['abs_error_lower'])
    assert ehi == F(row['abs_error_upper'])
    assert denominator*elo == F(row['primitive_error_lower'])
    assert denominator*ehi == F(row['primitive_error_upper'])
    if denominator*ehi < 1:
        small.append((n, r))
    else:
        assert denominator*elo > 1
assert small == [(1, 1), (1, 2), (2, 1)]
assert data[10][0][0] == 13657436403073

# Check the local denominator lemma through order four, at primes and powers.
def valuation(a, p):
    if not a:
        return None
    u, v, result = abs(a.numerator), a.denominator, 0
    while u % p == 0:
        u //= p
        result += 1
    while v % p == 0:
        v //= p
        result -= 1
    return result


local_checks = 0
for p in (2, 3, 5, 7):
    for e in (1, 2, 3):
        for offset in range(4):
            A = [F(p**e), F(1+p*offset)]+[F((j+offset)**2) for j in range(2, R+1)]
            b = [F(1+p*(offset+1))]+[F(j-offset) for j in range(1, R+1)]
            q = quotient(b, A)
            for r in range(R+1):
                assert valuation(q[r], p) == -(r+1)*e
                local_checks += 1

print('PASS: 14 independently constructed Newton-quotient jets; '
      '24 recurrence jet identities and 12 exact Casoratians through order four; '
      'all 20 saved rational/error/primitive certificates; exactly three '
      f'primitive errors below one; {local_checks} local valuation checks')
