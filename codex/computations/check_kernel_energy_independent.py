#!/usr/bin/env python3
"""Independent standard-library arithmetic and operator checks for the short energy proof."""
from fractions import Fraction as F
from math import comb, lcm

def harmonic(n,r): return sum((F(1,k**r) for k in range(1,n+1)),F(0))
def C(n,k): return comb(n,k)*comb(n+k,k)
def inner(a,b): return sum((u*v/F(i+j+1) for i,u in enumerate(a) for j,v in enumerate(b)),F(0))
def monomial_integral(i,j):
    # Pair (coefficient of zeta3, rational constant).
    if i==j: return F(2),-2*harmonic(i,3)
    return F(0),(harmonic(j,2)-harmonic(i,2))/F(j-i)

for n in range(1,26):
    L=[F((-1)**j*C(n,j)) for j in range(n+1)]
    residues=[F((-1)**(n-j)*C(n,j)) for j in range(n+1)]
    a=[c*c for c in residues]
    half_b=[c*sum((d/F(j-k) for j,d in enumerate(residues) if j!=k),F(0)) for k,c in enumerate(residues)]
    A=sum(a); B=sum(a[k]*harmonic(k,3)+half_b[k]*harmonic(k,2) for k in range(n+1))
    d=lcm(*range(1,n+1))
    assert sum(half_b)==0
    assert all((d*b).denominator==1 for b in half_b)
    assert (d**3*B).denominator==1
    # Recover quotient by triangular linear algebra on ordinary monomials.
    rem=L[:]; q=[F(0)]*n
    for j in reversed(range(n)):
        q[j]=rem[j+1]/(j+1)**2
        rem[j+1]-=q[j]*(j+1)**2
        rem[j]-=q[j]*(n*(n+1)-j*(j+1))
    assert all(v==0 for v in rem[1:]) and rem[0]==A
    assert inner(L,q)==0 and inner(L,L)==F(1,2*n+1)
    if n<=14:
        coeff=constant=F(0)
        for i,u in enumerate(L):
            for j,v in enumerate(L):
                c,b=monomial_integral(i,j);coeff+=u*v*c;constant+=u*v*b
        assert coeff==2*A and constant==-2*B
assert 64*3**9<F(169,5)**4
print('PASS: 25 exact residue/arithmetic and independent triangular-operator constructions; 14 independently evaluated logarithmic double integrals; exact exponential-rate comparison')
