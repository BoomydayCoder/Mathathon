#!/usr/bin/env python3
"""Exact checks and explicitly numerical evidence for multiplicity-lifted forms."""
from fractions import Fraction as F
from math import comb, factorial, isqrt
import mpmath as mp


def coefficients(n,m):
    return [(-1)**k*comb(m*n,k)*comb(n+k,k)**m for k in range(m*n+1)]


def form(n,m):
    a = coefficients(n,m)
    s = 2*m+1
    Hlow = [F(0)]
    Hhigh = [F(0)]
    for k in range(1,m*n+1):
        Hlow.append(Hlow[-1]+F(1,k**(s-1)))
        Hhigh.append(Hhigh[-1]+F(1,k**s))
    A = (s-1)*sum(v*v for v in a)
    B = sum(((s-1)*a[k]**2*Hhigh[k] for k in range(m*n+1)), F(0))
    B -= 2*sum((a[k]*a[l]*(Hlow[l]-Hlow[k])/F(l-k)
                for k in range(m*n+1) for l in range(k+1,m*n+1)),F(0))
    return A,B


def vp(a,p):
    result = 0
    while a and a % p == 0:
        result += 1
        a //= p
    return result


def residue(r,p):
    assert r.denominator % p
    return r.numerator*pow(r.denominator,-1,p) % p


def is_prime(p):
    return p >= 2 and all(p % d for d in range(2,isqrt(p)+1))


for m in range(1,5):
    for n in range(1,7):
        a = coefficients(n,m)
        for j in range(1,n+1):
            for r in range(m):
                derivative = sum((F((-1)**r*factorial(r)*a[k],(j+k)**(r+1))
                                  for k in range(m*n+1)),F(0))
                assert derivative == 0, (n,m,j,r)
        # Residue identity also checked at a non-root, non-pole point.
        t = m*n+2
        product = F((-1)**(m*n)*factorial(m*n),factorial(n)**m)
        for j in range(1,n+1): product *= (t-j)**m
        for k in range(m*n+1): product /= t+k
        assert sum((F(a[k],t+k) for k in range(m*n+1)),F(0)) == product

assert form(1,2) == (584,F(9689,16))
checks = 0
for m in [2,3,4]:
    for n in range(2,13):
        a = coefficients(n,m)
        A,B = form(n,m)
        s = 2*m+1
        for p in range(max(3,n+1,(m*n)//2+1),m*n+1):
            if not is_prime(p):continue
            r = m*n-p
            expected = (s+1)*sum(a[u]**2 for u in range(r+1)) % p
            assert residue(B*p**s,p) == expected
            checks += 1

for p in range(5,102):
    if not is_prime(p):continue
    n = (p+1)//2
    A,B = form(n,2)
    assert residue(B*p**5,p) == residue(F(291,8),p)
    if p != 97:
        assert vp(B.numerator,p)-vp(B.denominator,p) == -5
print('PASS exact Mellin zeros, rational transform, and prime-band congruences:', checks)

for n in [1,2,3,5,10,20,30]:
    A,B = form(n,2)
    precision = 15*n+80
    with mp.workdps(precision):
        first = A*mp.zeta(5)-mp.mpf(B.numerator)/B.denominator
    with mp.workdps(precision+60):
        second = A*mp.zeta(5)-mp.mpf(B.numerator)/B.denominator
        assert abs((first-second)/second) < mp.mpf('1e-35')
        an = 1/(mp.mpf(n)**3*comb(3*n,n))
        bound = mp.pi**mp.mpf('1.5')*an**2/6*(
            27*n**3*mp.gamma(2*n-mp.mpf('.5'))/mp.gamma(2*n+1)
            +9*n/mp.pi**2*mp.gamma(2*n+mp.mpf('.5'))/mp.gamma(2*n+1))
        assert abs(second) < bound
        print(n,mp.nstr(second,18),'logabs/n=',mp.nstr(mp.log(abs(second))/n,12))
print('Numerical values are evidence; exponential decay and nonvanishing are proved in the note.')
