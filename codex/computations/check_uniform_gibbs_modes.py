#!/usr/bin/env python3
"""Exact checks of the uniform-sampler modes and its differential inverse."""
from fractions import Fraction as F
from math import comb
import sympy as s

x=s.symbols('x')
for k in range(1,41):
    ln=[F((-1)**j*comb(k,j)*comb(k+j,j)) for j in range(k+1)]
    prev=[F((-1)**j*comb(k-1,j)*comb(k-1+j,j)) for j in range(k)]
    q=[];carry=F(0)
    for j in range(k):
        carry+=ln[j]+prev[j]
        q.append(carry)
    assert ln[k]==-q[-1]
    pq=[sum((q[j]/(j+1)**2 for j in range(r,k)),F(0)) for r in range(k)]
    assert pq==[c/k**2 for c in q]
for j in range(31):
    pf=sum(x**r for r in range(j+1))/s.Integer((j+1)**2)
    jj=-s.diff(x*(1-x)**2*s.diff(pf,x),x)/(1-x)
    assert s.cancel(jj+pf-x**j)==0
sig=s.symbols('sigma');kk=s.Function('K')(x)
hh=(1-x)*kk
ode=x*(1-x)*s.diff(kk,x,2)+(1-3*x)*s.diff(kk,x)-(1-sig)*kk
assert s.expand(s.diff(x*s.diff(hh,x),x)+sig*kk-ode)==0
print('PASS:40 exact uniform-sampler eigenmodes,31 polynomial inverse identities,')
print('     and the exact scalar ODE-to-Volterra-resolvent identity.')
