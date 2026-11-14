#!/usr/bin/env python3
"""Adversarial exact normalizer tests from direct node products.

Includes coefficient denominators with primes in (K,2K] and above 2K.
No floating point, harmonic diagnostics, or balanced-system solve is needed.
"""
from fractions import Fraction as Q
from math import factorial, gcd, lcm, prod


def prime(p):
    return p >= 2 and all(p % d for d in range(2, p) if d * d <= p)


def next_prime(p):
    while not prime(p):
        p += 1
    return p


def evaluate(coeff, t):
    out = 0
    for a in reversed(coeff):
        out = out * t + a
    return out


def affine(m, K, numerator):
    harmonic, a, b = Q(0), Q(0), Q(0)
    for j in range(1, K + 1):
        harmonic += Q(1, j ** (2 * m + 1))
        # Use the actual other-node product, not the factorial residue formula.
        den = prod(k*k - j*j for k in range(1, K + 1) if k != j)
        c = Q(numerator(-j*j), den)
        beta = -j**(2*m)*harmonic - Q(1, 2*m) + Q(1, 2*j)
        a += c*j**(2*m)
        b += c*beta
    return a, b


def valuation(n, p):
    n = abs(n)
    assert n
    result = 0
    while n % p == 0:
        n //= p
        result += 1
    return result


cases = inherited = 0
for m in range(1, 4):
    for d in range(4):
        K = 2*d + m + 1
        band_prime = next_prime(K + 1)
        assert band_prime <= 2*K
        denominators = [1] if d == 0 else [
            1, band_prime, band_prime**2, next_prime(2*K + 1)
        ]
        for T in denominators:
            P = [1] if d == 0 else [1] + [0]*(d-1) + [T]
            assert gcd(*P) == 1
            a, b = affine(m, K, lambda x: evaluate(P, x)**2)
            assert a == (-1)**m*T*T
            q = (-1)**(m+1)*b/(T*T)
            values = [evaluate(P, -j*j) for j in range(2*d+1)]
            G = gcd(*values)
            assert factorial(2*d) % G == 0
            diffs, leading = values, []
            while diffs:
                leading.append(diffs[0])
                diffs = [v-u for u, v in zip(diffs, diffs[1:])]
            assert gcd(*leading) == G
            assert leading[-1] == (-1)**d*T*factorial(2*d)
            L = lcm(*range(1, K+1))
            grid = Q(T*T*m*factorial(2*K)*L**(2*m+1), G*G)
            assert grid.denominator == 1 and (grid*q).denominator == 1
            band = prod(p for p in range(K+1, 2*K+1) if prime(p))
            refined = grid/band
            assert refined.denominator == 1 and (refined*q).denominator == 1
            for p in range(K+1, 2*K+1):
                if prime(p):
                    assert valuation(refined.numerator, p) == 2*valuation(T, p)
                    assert valuation(q.denominator, p) <= 2*valuation(T, p)
            remaining = q.denominator
            for p in range(2, max(K, T)+1):
                if prime(p) and (p <= K or T % p == 0):
                    while remaining % p == 0:
                        remaining //= p
            assert remaining == 1
            if T > 1 and gcd(q.denominator, T) > 1:
                inherited += 1
            cases += 1

# Properness alone controls affine coefficients, not division by their slope.
a, b = affine(1, 2, lambda x: x+2)
assert (a, b) == (3, Q(-7, 2))
assert -b/a == Q(7, 6)

# The polynomial quotient in an improper numerator can restore an upper prime.
a, b = affine(1, 2, lambda x: -85*x-84)
quotient_moment = Q(1, 12)-5*Q(1, 12)+21*Q(1, 4)
assert (a, b + quotient_moment) == (-341, Q(4919, 12))

# One saved balanced example is reconstructed without its moment-matrix data.
a, b = affine(1, 4, lambda x: (3*x-1)**2)
assert a == -9 and b/9 == Q(84127, 69984)
assert b/9 - Q(577, 480) == Q(1, 174960)
assert inherited > 0
print(f"PASS: {cases} arbitrary rational trials, including {inherited} inherited "
      "denominator-prime cases; affine-root and improper counterexamples; exact balanced example")
