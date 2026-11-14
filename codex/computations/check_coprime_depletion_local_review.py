"""Independent exact determinant checks using cyclotomic cumulants.

Polynomial moments are generated from exp(2z)-exp(z)+1, not Bernoulli
numbers. Pole constants are generated from the divisibility-weight sum,
not from the four-Hurwitz-tail formula in the construction script.
"""
from fractions import Fraction as F
from math import comb
import sympy as S

t, X = S.symbols('t X')


def valuation(q, p):
    n, d = map(int, S.fraction(q))
    if not n:
        return float('inf')
    def count(v):
        v = abs(v)
        result = 0
        while v % p == 0:
            v //= p
            result += 1
        return result
    return count(n)-count(d)


def sigma_polynomial(s, e):
    k = 2*e+2
    exponential = [1]+[2**j-1 for j in range(1, k+1)]
    cumulants = [0]
    for n in range(1, k+1):
        cumulants.append(exponential[n]-sum(comb(n-1, j-1)*cumulants[j]
                                            *exponential[n-j] for j in range(1, n)))
    return S.Rational((-1)**e*cumulants[k]*comb(2*e+s, s-1), 6**k)


def sigma_pole(s, r):
    m = s-1
    weighted = sum((F((1-2*int(n % 2 == 0))*(1-3*int(n % 3 == 0)), n**s)
                    for n in range(1, 6*r+1)), F(0))
    constant = -(6*r)**m*weighted+F(1, 6*r)
    leading = (2**m-1)*(3**m-1)*r**m
    return leading*X+S.Rational(constant.numerator, constant.denominator)


def determinant(s, K, N):
    h = K-N
    numerator = S.Poly(S.prod(t+j*j for j in range(1, N+1)), t)**6
    denominator = S.Poly(S.prod(t+j*j for j in range(1, K+1)), t)
    moments = []
    for d in range(2*h-1):
        quotient = S.quo(numerator*S.Poly(t**d, t), denominator)
        value = sum(c*sigma_polynomial(s, e) for (e,), c in quotient.terms() if c)
        for r in range(N+1, K+1):
            residue = S.Rational(numerator.eval(-r*r)*(-r*r)**d,
                                 S.prod(j*j-r*r for j in range(1, K+1) if j != r))
            value += residue*sigma_pole(s, r)
        moments.append(S.expand(value))
    result = S.Poly(S.Matrix(h, h, lambda i, j: moments[i+j]).det(method='domain-ge'), X)
    expected_leading = ((-1)**(h*(h-1)//2)*((2**(s-1)-1)*(3**(s-1)-1))**h
                       *S.prod(r**(s-1)*S.prod(j*j-r*r for j in range(1, N+1))**5
                               for r in range(N+1, K+1)))
    assert result.degree() == h and result.LC() == expected_leading
    return result


if __name__ == '__main__':
    cases = ((11, 1, 0), (3, 6, 0), (5, 6, 0), (7, 6, 0),
             (5, 7, 1), (7, 8, 3), (9, 6, 5), (9, 4, 0))
    checked = 0
    for s, K, N in cases:
        P = determinant(s, K, N)
        h = K-N
        fixed = (2**(s-1)-1)*(3**(s-1)-1)*(2**s-1)
        M, primitive = P.monic().clear_denoms(convert=True)
        assert primitive.primitive()[0] == 1
        count = 0
        for p in S.primerange(max(3, 2*K)+1, 6*K+1):
            if fixed % p == 0:
                continue
            b = K-max(N, p//6)
            assert valuation(P.LC(), p) == 0
            assert min(valuation(c, p) for c in P.all_coeffs()) == -s*b
            assert valuation(P.nth(h-b), p) == -s*b
            assert valuation(M, p) == s*b
            checked += 1
            count += 1
        print(dict(s=s, K=K, N=N, degree=h, admissible_primes=count), flush=True)
    # A fixed-prime exception is a real cancellation, not just a presentation choice.
    exceptional = determinant(3, 3, 0)
    assert min(valuation(c, 7) for c in exceptional.all_coeffs()) == -5
    assert -3*(3-7//6) == -6
    print(f'PASS: {checked} new determinant-prime checks, plus the p=7,s=3 exception.')
