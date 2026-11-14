#!/usr/bin/env python3
"""Construct Wilson minimizers from moments, without a hypergeometric formula.

Checks the proved moving-pole family against Apéry recurrence values. This
does not test or assume the conjectural D_(2h) norm-duality identity.
"""
import sympy as sp

t, x, y = sp.symbols("t x y")


def moment(k):
    return (-1) ** k * sp.Rational(2 * k + 3, 2) * sp.bernoulli(2 * k + 2)


def mu(poly):
    return sum(c * moment(e[0]) for e, c in sp.Poly(sp.expand(poly), t).terms())


def functional(poly):
    return sum(
        c * e[0] * sp.bernoulli(e[0] - 1, 0)
        for e, c in sp.Poly(sp.expand(poly), x).terms()
        if e[0]
    )


def bound(poly, a):
    r = mu(sp.cancel(t * (1 - poly / poly.subs(t, -a * a)) / (t + a * a)))
    return (
        sp.harmonic(a, 3)
        + sp.Rational(1, 2 * a * a)
        - sp.Rational(1, 2 * a**3)
        + sp.Rational(1, 4 * a**4)
        - r / a**4
    )


A, B = [sp.Integer(1), sp.Integer(5)], [sp.Integer(0), sp.Integer(6)]
for n in range(1, 8):
    beta = (2 * n + 1) * (17 * n * n + 17 * n + 5)
    for seq in (A, B):
        seq.append((beta * seq[-1] - n**3 * seq[-2]) / (n + 1) ** 3)

for n in range(8):
    # Orthogonality for dν=t dμ, with monic leading coefficient fixed.
    if n:
        gram = sp.Matrix(n, n, lambda i, j: moment(i + j + 1))
        rhs = sp.Matrix([-moment(n + i + 1) for i in range(n)])
        coeff = gram.inv() * rhs
        poly = t**n + sum(coeff[j] * t**j for j in range(n))
    else:
        poly = sp.Integer(1)
    assert all(mu(t ** (j + 1) * poly) == 0 for j in range(n))
    h = n + 1
    pi = sp.expand(sum(
        sp.binomial(h, k) * sp.binomial(h + k, k)
        * sp.prod((y - j) * (y + j + 1) for j in range(k)) / sp.factorial(k)**2
        for k in range(h + 1)
    ))
    delta = sp.Poly(sp.cancel((pi - pi.subs(y, y - 1)) / (2 * y)), y)
    transformed = sp.Poly(sum(c * (-t) ** (e[0] // 2) for e, c in delta.terms()), t)
    assert all(e[0] % 2 == 0 for e, c in delta.terms())
    assert sp.expand(poly - transformed.monic().as_expr()) == 0
    assert bound(poly, h) == (5 * B[h] - B[h - 1]) / (5 * A[h] - A[h - 1])

# Include nonzero constant terms: the correction functional is then essential.
for degree in range(1, 7):
    pi = sp.expand(7 + sum((-1)**j * (j + 2) * (y * (y + 1))**j
                           for j in range(1, degree + 1)))
    delta = sp.expand(pi - pi.subs(y, y - 1))
    d = sp.Poly(sp.cancel(delta / (2 * y)), y)
    poly = sum(c * (-t) ** (e[0] // 2) for e, c in d.terms())
    theta = functional(sp.cancel((pi.subs(y, x) - pi) / (x - y)))
    correction = functional(x * pi.subs(y, x))
    for a in (1, 2, 5):
        den = delta.subs(y, a)
        assert den != 0
        associated = (
            sp.harmonic(a, 3)
            + (pi.subs(y, a - 1) / a**3
               + (theta.subs(y, a) - theta.subs(y, a - 1)) / 2) / den
        )
        assert sp.cancel(associated - bound(poly, a) - correction / (a**3 * den)) == 0

print("PASS: eight moment-constructed minimizers and adjacent-Apéry bounds; "
      "18 mixed-basis correction identities including constant terms")
