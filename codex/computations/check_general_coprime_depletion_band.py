"""Exact step-value, band-integral and determinant checks for coprime depletion."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, gcd
import sympy as S

t, z, X = S.symbols('t z X')


def vp(q, p):
    num, den = map(int, S.fraction(q))
    if not num:
        return float('inf')
    def count(n):
        n = abs(n)
        out = 0
        while n % p == 0:
            n //= p
            out += 1
        return out
    return count(num)-count(den)


def harmonic(k, s):
    return sum((F(1, j**s) for j in range(1, k+1)), F(0))


@lru_cache(None)
def steps(a, b, s):
    q, m = a*b, s-1
    return tuple(q**m*harmonic(k, s)-a**m*harmonic(k//b, s)
                 -b**m*harmonic(k//a, s) for k in range(1, (q-1)//2+1))


@lru_cache(None)
def cumulant(a, b, k):
    q = a*b
    polynomial = S.Poly(S.cancel((z**q-1)*(z-1)/((z**a-1)*(z**b-1))), z)
    assert polynomial.domain == S.ZZ and polynomial.eval(1) == 1
    exponential = [sum(int(c)*e**n for (e,), c in polynomial.terms())
                   for n in range(k+1)]
    values = [0]
    for n in range(1, k+1):
        values.append(exponential[n]-sum(comb(n-1, j-1)*values[j]*exponential[n-j]
                                         for j in range(1, n)))
    return values[k]


def moment(a, b, s, e):
    k = 2*e+2
    return S.Rational((-1)**e*cumulant(a, b, k)*comb(2*e+s, s-1), (a*b)**k)


def pole(a, b, s, r):
    q, m = a*b, s-1
    weighted = sum((F((1-a*int(k % a == 0))*(1-b*int(k % b == 0)), k**s)
                    for k in range(1, q*r+1)), F(0))
    constant = -(q*r)**m*weighted+F((a-1)*(b-1), 2*q*r)
    return ((a**m-1)*(b**m-1)*r**m*X
            +S.Rational(constant.numerator, constant.denominator))


def determinant(a, b, s, K, N):
    h = K-N
    numerator = S.Poly(S.prod(t+j*j for j in range(1, N+1)), t)**6
    denominator = S.Poly(S.prod(t+j*j for j in range(1, K+1)), t)
    values = []
    for d in range(2*h-1):
        quotient = S.quo(numerator*S.Poly(t**d, t), denominator)
        value = sum(c*moment(a, b, s, e) for (e,), c in quotient.terms() if c)
        for r in range(N+1, K+1):
            residue = S.Rational(numerator.eval(-r*r)*(-r*r)**d,
                                 S.prod(j*j-r*r for j in range(1, K+1) if j != r))
            value += residue*pole(a, b, s, r)
        values.append(S.expand(value))
    return S.Poly(S.Matrix(h, h, lambda i, j: values[i+j]).det(method='domain-ge'), X)


def B(q, alpha):
    if alpha <= F(2, q):
        return F((q-2)**2, 2*q)
    return (1-alpha)*(q*(1+alpha)/2-2)


if __name__ == '__main__':
    step_checks = 0
    integral_checks = 0
    for a in range(2, 9):
        for b in range(a+1, 9):
            if gcd(a, b) != 1:
                continue
            q = a*b
            for s in (3, 5, 7):
                for value in steps(a, b, s):
                    assert value > F(79, 144)*q**(s-1)
                    step_checks += 1
            for alpha in (F(0), F(3, 40), F(1, 3), F(1, 2), F(9, 10), F(2, q)):
                cuts = sorted({F(2), F(q)} | ({q*alpha} if 2 < q*alpha < q else set()))
                direct = sum((v-u)*((1-max(alpha, u/q))+(1-max(alpha, v/q)))/2
                             for u, v in zip(cuts, cuts[1:]))
                assert direct == B(q, alpha) and B(q, alpha) >= B(6, alpha)
                integral_checks += 1
    cases = ((2, 5, 3, 3, 0), (3, 4, 5, 3, 1), (3, 5, 3, 4, 1),
             (4, 5, 5, 3, 0), (4, 7, 3, 2, 1), (5, 6, 7, 2, 0))
    determinant_checks = 0
    for a, b, s, K, N in cases:
        q, m, h = a*b, s-1, K-N
        C = (a**m-1)*(b**m-1)
        P = determinant(a, b, s, K, N)
        leading = ((-1)**(h*(h-1)//2)*C**h
                   *S.prod(r**m*S.prod(j*j-r*r for j in range(1, N+1))**5
                           for r in range(N+1, K+1)))
        assert P.degree() == h and P.LC() == leading
        M, primitive = P.monic().clear_denoms(convert=True)
        assert primitive.primitive()[0] == 1
        count = 0
        for p in S.primerange(2*K+1, q*K+1):
            if 2*q*C % p == 0 or any(vp(S.Rational(v.numerator, v.denominator), p) != 0
                                     for v in steps(a, b, s)):
                continue
            bad = K-max(N, p//q)
            assert min(vp(c, p) for c in P.all_coeffs()) == -s*bad
            assert vp(P.nth(h-bad), p) == -s*bad
            assert vp(M, p) == s*bad
            count += 1
            determinant_checks += 1
        print(dict(a=a, b=b, s=s, K=K, N=N, checked_primes=count), flush=True)
    print(f'PASS: {step_checks} positive step values; {integral_checks} exact band integrals; '
          f'{determinant_checks} determinant-prime checks.')
