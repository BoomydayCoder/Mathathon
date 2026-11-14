#!/usr/bin/env python3
"""Independent finite checks from the unreduced rational products.

The general proof is in apery-mixed-numerator-independent-review.md.
These checks include asymmetric boundary cases and both summation endpoints.
"""
import sympy as sp

t = sp.symbols("t")


def product(n):
    return sp.prod(t - j for j in range(1, n + 1)) / sp.prod(
        t + j for j in range(n + 1)
    )


def mixed(n, m):
    return (-1) ** (n + m) * product(n) * product(m)


for n in range(5):
    for m in range(5):
        q = sp.cancel((t - m - 1) * (t + m) * mixed(n, m))
        lhs = (n + 1) ** 3 * mixed(n + 1, m)
        lhs -= (2 * n + 1) * (n * (n + 1) + 2 * m * (m + 1) + 1) * mixed(n, m)
        if n:
            lhs += n**3 * mixed(n - 1, m)
        assert sp.cancel(lhs + 2 * (2 * n + 1) * (q.subs(t, t + 1) - q)) == 0
        assert sp.diff(q, t).subs(t, 1) == (1 if n == 0 else 0)
        assert sp.limit(sp.diff(q, t), t, sp.oo) == 0

for n in range(1, 8):
    certificate = sp.cancel(4 * (t + n) ** 2 * product(n) ** 2 / n**2)
    lhs = product(n) ** 2 + product(n - 1) ** 2 + 6 * product(n) * product(n - 1)
    assert sp.cancel(lhs - certificate.subs(t, t + 1) + certificate) == 0
    assert sp.diff(certificate, t).subs(t, 1) == 0
    assert sp.limit(sp.diff(certificate, t), t, sp.oo) == 0

print("PASS: 25 asymmetric mixed products, seven adjacent products, all derivative endpoints")
