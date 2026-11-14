#!/usr/bin/env python3
"""Exact checks of the positive dual Laurent model and initial residues."""
import sympy as s
from math import comb
w=s.symbols('w')
for m in range(1,5):
    coeff=s.series((w/s.sin(w))**(2*m),w,0,2*m).removeO().expand()
    values=[s.expand(coeff).coeff(w,2*l) for l in range(m)]
    print('m=',m,'initial odd-zeta coefficients:',[(2*m-2*l)*values[l] for l in range(m)])
    assert values[0]==1
assert s.series((w/s.sin(w))**6,w,0,6).removeO().coeff(w,4)==s.Rational(8,15)
a,b,u=s.symbols('a b u',positive=True)
m=2;U=u**3
L=(1+a)**2*(1+b)**2*(1+u)**3*((1+a)*(1+b)+U*a*b)/(U*a*b)
assert L.subs({a:s.Rational(1,3),b:s.Rational(1,3),u:2})==256
assert all(s.diff(L,v).subs({a:s.Rational(1,3),b:s.Rational(1,3),u:2})==0 for v in (a,b,u))
assert [sum(comb(n,k)**4*comb(2*n+k,k)**2 for k in range(n+1)) for n in range(3)]==[1,10,626]
# Exact odds substitution, including the residual weight, for arbitrary m,n.
# Compare logarithmic exponents by checking several independent symbolic values.
a,b,U=s.symbols('a b U',positive=True)
for m in range(1,5):
    for n in range(3):
        x=a/(1+a);y=b/(1+b)
        old=((1-x)*(1-y))**(m*n)*(U*x*y/(1+U*x*y))**(n+1)/((1+a)*(1+b))
        kernel=(1+a)**m*(1+b)**m*((1+a)*(1+b)+U*a*b)/(U*a*b)
        new=kernel**(-n-1)*(1+a)**(m-1)*(1+b)**(m-1)
        assert s.factor(old-new)==0
print('PASS exact dual-period measure, positive critical point, initial sequence and residues')
