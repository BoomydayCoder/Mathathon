#!/usr/bin/env python3
"""Focused exact checks supplementing the adaptive-boundary proof audit."""
from fractions import Fraction as F
from math import comb,factorial,lcm
import sympy as S

for n in range(5,201):
    g2=comb(n,2)*comb(n+2,2)
    assert comb(2*n,n)>=g2>=F(n**4,4)
    ratios=[F((n+j+1)*(n-j),(j+1)**2) for j in range(2,n)]
    assert all(a>b for a,b in zip(ratios,ratios[1:]))

for n in range(1,31):
    Dn=lcm(*range(1,n+1));D2n=lcm(*range(1,2*n+1))
    K=Dn**2*lcm(*range(1,2*n+2))**4
    harmonic=[F(0)]
    for k in range(1,2*n+1):harmonic.append(harmonic[-1]+F(1,k))
    D=E=G=F(0)
    for j in range(1,n+1):
        u=F(factorial(j-1)**2*factorial(n-j),factorial(n+j))
        assert (Dn*D2n*u).denominator==1
        assert abs(u)==F(1,j*j*comb(n,j)*comb(n+j,j))
        delta=harmonic[n-j]+harmonic[n+j]-2*harmonic[j-1]
        assert abs(delta)<=4*harmonic[2*n]
        D+=u*u*delta
        E+=u*u*delta/F(j+n+1)
        G+=u*u/F((j+n+1)**2)
    W=2*E+G
    assert (K*D).denominator==(K*W).denominator==1
    assert D>0 and W>0  # finite supplement; general proof is eventual
    a=2*D/W-(n+1)
    d=D/F(n+1)/W
    assert 0<a<1 and 0<d<1

n,a,y=S.symbols('n a y',real=True)
c=n+S.Rational(1,2);b=n+1;z=c+S.I*y
assert S.cancel((z-a)/(z+b)-S.Rational(1,2)
                -(-S.Rational(1,2)-2*a+S.I*y)/(2*(c+b+S.I*y)))==0
print('PASS uniform endpoint comparisons, exact root weights, arithmetic clearing and contour multiplier')
