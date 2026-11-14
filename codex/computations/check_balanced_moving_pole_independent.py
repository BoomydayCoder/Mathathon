#!/usr/bin/env python3
"""Independent arithmetic check: ordinary quotient moments vs reflected Bernoulli formula."""
from fractions import Fraction as F
from math import comb,factorial,lcm
from functools import lru_cache
import json
import sympy as S
from pathlib import Path
rows=json.loads((Path(__file__).parent/'balanced_rodrigues_integrality_checks.json').read_text())
@lru_cache(None)
def bern(k):
    b=S.bernoulli(k) if k!=1 else -S.Rational(1,2)
    return F(int(S.numer(b)),int(S.denom(b)))
def moment(m,k):return (-1)**k*F(comb(2*k+2*m+1,2*m),2*k+2)*bern(2*k+2)
def evaluate(c,x):return sum(v*x**k for k,v in enumerate(c))
def quotient(c,a):
    # Quotient of c(x)-c(a) by x-a, via descending monomials.
    out=[F(0)]*(len(c)-1);carry=F(0)
    for k in range(len(c)-1,0,-1):
        carry=c[k]+a*carry;out[k-1]=carry
    assert c[0]+a*carry==evaluate(c,a)
    return out
count=0
for row in rows:
    m,n=row['m'],row['n'];d=m*n
    q=[F(v) for v in row['integer_polynomial_coefficients_low_to_high']]
    normalized=[v/q[0] for v in q];H=lcm(*range(1,n+1))**(2*m-2)
    f=[F(0)]*(2*d+1)
    for k,v in enumerate(normalized):f[2*k]=H*(-1)**k*v
    for a in [1,n+1,2*n+3]:
        fa=evaluate(f,a);assert fa.denominator==1
        g=quotient(normalized,-a*a)
        mu_g=sum(c*moment(m,k) for k,c in enumerate(g))
        P=[F(0)]*(2*m)+quotient(f,a)
        derivative=2*m-1
        bb=sum(c*F(factorial(k),factorial(k-derivative))*bern(k-derivative)
               for k,c in enumerate(P) if k>=derivative)
        assert H*mu_g==-bb/factorial(2*m)-(fa-H)/F(2*a)
        harmonic=sum((F(1,j**(2*m+1)) for j in range(1,a+1)),F(0))
        rational=-a**(2*m)*fa*harmonic-fa/F(2*m)+H/F(2*a)-bb/factorial(2*m)
        direct=H*mu_g+fa*(-a**(2*m)*harmonic-F(1,2*m)+F(1,2*a))
        assert rational==direct
        J=max(a,2*d+2*m-1);M=factorial(2*m)*lcm(*range(1,J+1))**(2*m+1)
        assert (M*rational).denominator==1
        assert (M*a**(2*m)*fa).denominator==1
        count+=1
print(f'PASS: {count} exact independent quotient-moment/reflection comparisons and full associated affine integer normalizers')
