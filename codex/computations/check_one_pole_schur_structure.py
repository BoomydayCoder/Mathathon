#!/usr/bin/env python3
"""Exact checks for the all-degree proof in one-pole-schur-structure.md."""
import sympy as S

t = S.symbols('t')


def moment(e, s=3):
    return ((-1)**e*S.bernoulli(2*e+2)*S.factorial(2*e+s)
            /(S.factorial(s-1)*S.factorial(2*e+2)))


def integral(p, shift=0):
    return sum(c*moment(mon[0]+shift)
               for mon,c in S.Poly(S.expand(p), t).terms())


for m in range(1, 16):
    n = m-1
    A = S.Matrix(n,n,lambda i,j:moment(i+j+1)+moment(i+j))
    B = S.Matrix(n,n,lambda i,j:moment(i+j+2)+moment(i+j+1))
    b = S.Matrix([moment(i) for i in range(n)])
    c = S.Matrix([moment(i+1) for i in range(n)])
    L = 1+(b.T*A.inv()*b)[0] if n else S.Integer(1)
    U = S.Rational(5,4)-(c.T*B.inv()*c)[0] if n else S.Rational(5,4)
    H = sum(S.Rational(1,j**3) for j in range(1,m+1))
    assert L == H+S.Rational(m-1,2*m**3)
    assert U == H+S.Rational(1,2*m*(m+1))
    assert U-L == S.Rational(1,2*m**3*(m+1))

for shift in [0,1]:
    for n in range(7):
        pref = S.rf(2,n)**(2+shift)*S.rf(1,n)**(1-shift)
        q = S.expand(sum(
            S.rf(-n,k)*S.rf(n+2+shift,k)
            *S.prod(t+j*j for j in range(1,k+1))
            /(S.rf(2,k)**(2+shift)*S.rf(1,k)**(1-shift)*S.factorial(k))
            for k in range(n+1)))
        assert q.subs(t,-1) == 1
        W = S.expand(pref*q)
        norm = (S.factorial(n+1)**2*S.factorial(n)**4/S.Integer(4*(n+1))
                if shift == 0 else
                S.factorial(n)*S.factorial(n+1)**6
                /(2*(2*n+3)*S.factorial(n+2)))
        assert integral(W**2,shift) == norm
        assert all(integral(W*t**j,shift) == 0 for j in range(n))
        if n:
            p = S.cancel((1-q)/(t+1))
            assert S.denom(p) == 1 or S.denom(p).is_number
            R = integral(p,shift)
            m = n+1
            H = sum(S.Rational(1,j**3) for j in range(1,m+1))
            expected = (H+S.Rational(m-1,2*m**3)-1 if shift == 0 else
                        S.Rational(5,4)-H-S.Rational(1,2*m*(m+1)))
            assert R == expected

m = S.symbols('m', positive=True, integer=True)
assert S.factor(1/(m+1)**3+m/(2*(m+1)**3)-(m-1)/(2*m**3)
                -(2*m+1)/(2*m**3*(m+1)**3)) == 0
assert S.factor(1/(2*m*(m+1))-1/(m+1)**3-1/(2*(m+1)*(m+2))
                -1/(m*(m+1)**3*(m+2))) == 0
for e in range(15):
    assert moment(e,5) == S.Rational((2*e+4)*(2*e+5),12)*moment(e,3)
print('PASS Schur identities m<=15; Wilson orthogonality and norms n<=6; exact increments and moment transform')

# Exact least-denominator checks, not checks of a merely sufficient normalizer.
from fractions import Fraction
from math import isqrt


def primes_through(n):
    flags = bytearray(b'\x01')*(n+1)
    flags[:2] = b'\x00\x00'
    for p in range(2,isqrt(n)+1):
        if flags[p]:
            flags[p*p:n+1:p] = b'\x00'*((n-p*p)//p+1)
    return [p for p in range(2,n+1) if flags[p]]


def vp_integer(a,p):
    exponent = 0
    while a and a % p == 0:
        exponent += 1
        a //= p
    return exponent


harmonic = Fraction(0)
primes = primes_through(300)
valuation_checks = exceptional_upper_checks = 0
for m in range(1,301):
    harmonic += Fraction(1,m**3)
    lower = harmonic+Fraction(m-1,2*m**3)
    upper = harmonic+Fraction(1,2*m*(m+1))
    for p in primes:
        if p >= 3 and m < 2*p and p < m:
            for endpoint in [lower,upper]:
                valuation = (vp_integer(endpoint.numerator,p)
                             -vp_integer(endpoint.denominator,p))
                assert valuation == -3, (m,p,endpoint)
                assert vp_integer(endpoint.denominator,p) == 3
                valuation_checks += 1
            if m+1 == 2*p:
                exceptional_upper_checks += 1
    if m in [10,30,100,300]:
        # Guaranteed lower bound for the primitive cleared error, from (28).
        assert Fraction(lower.denominator,4*(m+1)**4) > 1
        assert Fraction(upper.denominator,4*(m+2)**4) > 1
print('PASS exact primitive denominator valuations:',valuation_checks,
      '; upper correction exceptional locations checked:',exceptional_upper_checks)
