#!/usr/bin/env python3
"""Exact polynomial checks and independent scalar integral corroboration."""
from fractions import Fraction as F
from math import comb
import mpmath as mp

def leg(n):return [F((-1)**j*comb(n,j)*comb(n+j,j)) for j in range(n+1)]
def P(c):return [sum((c[j]/(j+1)**2 for j in range(k,len(c))),F(0)) for k in range(len(c))]
def pair(c,d):return sum((x*y/F((i+j+1)*(i+j+2)) for i,x in enumerate(c) for j,y in enumerate(d)),F(0))
for k in range(1,41):
    a,b=leg(k),leg(k-1);q=[];carry=F(0)
    for j in range(k):
        carry+=a[j]+b[j];q.append(carry)
    assert P(q)==[v/k**2 for v in q]
    assert pair(q,q)==F(2,k)
for i in range(15):
    x=[F(0)]*i+[F(1)]
    for j in range(15):
        y=[F(0)]*j+[F(1)]
        assert pair(x,P(y))==F(1,(i+1)*(j+1)*(i+j+2))
        assert pair(x,P(y))==pair(P(x),y)
print('PASS 40 exact mode eigenvalues/norms and 225 independent weighted-Hilbert pairings')
mp.mp.dps=60
f=lambda t:(mp.zeta(2)-mp.polylog(2,1-t))**2/t if t else mp.mpf(0)
value=mp.quad(f,[0,mp.mpf('.1'),mp.mpf('.5'),1])
error=abs(value-2*mp.zeta(5))
assert error<mp.mpf('1e-55')
print('PASS direct dilogarithm positive integral for 2*zeta(5), error',mp.nstr(error,5))
