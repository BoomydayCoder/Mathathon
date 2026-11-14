#!/usr/bin/env python3
"""Independent diagonal numerator checks from ordinary rational moments."""
from fractions import Fraction as F
from math import comb, factorial, lcm
from functools import lru_cache
from pathlib import Path
import json
import sympy as S
@lru_cache(None)
def bern(k):
    z=S.bernoulli(k) if k!=1 else -S.Rational(1,2)
    return F(int(S.numer(z)),int(S.denom(z)))
def moment(m,k):
    return (-1)**k*F(comb(2*k+2*m+1,2*m),2*k+2)*bern(2*k+2)
def at(c,x):return sum(v*x**k for k,v in enumerate(c))
def quotient(c,x):
    q=[F(0)]*(len(c)-1);v=F(0)
    for k in range(len(c)-1,0,-1):
        v=c[k]+x*v;q[k-1]=v
    assert c[0]+x*v==at(c,x)
    return q
rows=json.loads(Path(__file__).with_name('balanced_rodrigues_integrality_checks.json').read_text())
for row in rows:
    m,n=row['m'],row['n'];a=n+1
    c=[F(x) for x in row['integer_polynomial_coefficients_low_to_high']]
    c=[x/c[0] for x in c]
    Ra=at(c,-a*a);assert Ra.denominator==1
    alpha=a**(2*m)*Ra
    q=quotient(c,-a*a)
    qmoment=sum(v*moment(m,k) for k,v in enumerate(q))
    hs=sum((F(1,j**(2*m+1)) for j in range(1,a+1)),F(0))
    b=alpha*hs+Ra*F(1,2*m)-Ra*F(1,2*a)-qmoment
    scaled=factorial(2*m)*b
    for prime,exponent in S.factorint(scaled.denominator).items():
        assert prime<=a
        if prime==a:assert exponent<=1
    L=2*m*n+2*m-1
    lam=1
    for prime in S.primerange(1,n+1):
        power=prime
        while power*prime<=L:power*=prime
        lam*=power
    U=factorial(2*m)*a*lcm(*range(1,n+1))**(2*m-2)*lam**(2*m+1)
    assert (U*b).denominator==1
    assert (U*alpha).denominator==1
print(f'PASS: {len(rows)} independently reconstructed diagonal coefficients, prime supports, and full normalizers')
