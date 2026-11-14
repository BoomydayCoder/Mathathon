#!/usr/bin/env python3
"""Independent coefficient-array Rodrigues check against Bernoulli moments."""
from fractions import Fraction as F
from math import comb, factorial
from functools import lru_cache
import sympy as S

def multiply(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):out[i+j]+=u*v
    return out

def raise_even(coeff,m,s):
    # coeff[k] multiplies x^(2k). Odd part of A(x) f(x+1/2), divided by x.
    shifted=[F(0)]*(2*len(coeff)-1)
    for k,c in enumerate(coeff):
        for j in range(2*k+1):shifted[j]+=c*comb(2*k,j)*F(1,2)**(2*k-j)
    N=2*m+1;a=F(s,2)
    power=[F(comb(N,j))*(a+1)**(N-j) for j in range(N+1)]
    product=multiply(multiply([a,F(1)],power),shifted)
    return product[1::2]

@lru_cache(None)
def moment(s,e):
    b=S.bernoulli(2*e+2)
    return F((-1)**e*comb(2*e+s,s-1)*int(S.numer(b)),(2*e+2)*int(S.denom(b)))

systems=conditions=0
for m in range(1,7):
    for n in range(1,7):
        coeff=[F(1)]
        for s in reversed(range(n)):coeff=raise_even(coeff,m,s)
        d=m*n;lead=F(factorial(d+m+n)//factorial(d+m))
        assert len(coeff)==d+1 and coeff[-1]==lead
        p=[(-1)**(d+k)*c/lead for k,c in enumerate(coeff)]
        assert p[-1]==1
        assert p[0]==F((-1)**d*factorial(n)**(2*m),comb((m+1)*n+m,n))
        assert sum(c*(-1)**k for k,c in enumerate(p))/p[0]==(n+1)**(2*m)
        for j in range(1,m+1):
            for k in range(n):
                assert sum(c*moment(2*j+1,e+k+m-j) for e,c in enumerate(p))==0
                conditions+=1
        systems+=1
print(f'PASS: {systems} independent coefficient-array constructions; {conditions} exact Bernoulli orthogonality conditions; both all-index endpoint formulas checked')
