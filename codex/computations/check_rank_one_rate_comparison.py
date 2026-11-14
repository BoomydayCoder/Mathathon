#!/usr/bin/env python3
"""Exact algebraic certificates plus labeled numerical corroboration."""
from math import factorial
from fractions import Fraction as Q
from pathlib import Path
import json
import sympy as S
import mpmath as mp
z=S.symbols('z')
A={1:S.Poly(1,z)}
for d in range(2,21):
    prev=A[d-1].as_expr()
    A[d]=S.Poly(S.expand((1+(d-1)*z)*prev+z*(1-z)*S.diff(prev,z)),z)
pairs=0
for m in range(1,11):
    poly=A[2*m];degree=2*m-1
    co=[poly.nth(i) for i in range(degree+1)]
    assert co==co[::-1] and sum(co)==factorial(2*m)
    assert co[0]==1 and co[1]>=1 and all(v>=0 for v in co)
    for k in range(m):
        quotient,remainder=S.div(1+z-z**k-z**(degree-k),1-z,z)
        assert remainder==0 and all(v>=0 for v in S.Poly(quotient,z).all_coeffs())
        pairs+=1
print('PASS:',pairs,'exact Eulerian pair inequalities; normalization and lower bound through s=21')

count=0
for h in range(1,9):
 for r in range(0,9):
  for y2 in (Q(0),Q(1,4),Q(1),Q(7),Q(1000)):
   value=Q(1)
   for j in range(2*h+1,2*h+r+1):value*=(1+y2)/(j*j+y2)
   lower=Q(factorial(2*h),factorial(2*h+r))**2
   assert lower<=value<=1
   count+=1
print('PASS:',count,'exact rational measure-factor comparisons including h=1 and s=3')

mp.mp.dps=160
for m in range(1,8):
    r=m-1;a=2*mp.pi
    c=mp.mpf(2)/factorial(2*m)*(a*a/(1+a*a))**r
    C=(1+a*a)**r
    for y in [mp.mpf('1e-8'),mp.mpf('.001'),mp.mpf('.1'),mp.mpf(1),mp.mpf(10),mp.mpf(100)]:
        zz=mp.exp(-a*y)
        T=y/(-mp.expm1(-a*y))
        polynomial=sum(int(A[2*m].nth(i))*zz**i for i in range(2*m))
        ratio=2*a**(2*r)/factorial(2*m)*T**(2*r)*polynomial/(1+zz)
        eps=mp.mpf('1e-140')
        assert c*(1+y*y)**r*(1-eps)<=ratio<=C*(1+y*y)**r*(1+eps)
print('PASS numerical corroboration: density bounds over seven odd indices and six scales')

rows=json.loads(Path(__file__).with_name('rank_one_subspace_approximants.json').read_text())
errors={}
for row in rows:
    s,h=row['s'],row['h']
    root=mp.mpf(row['numerator'])/mp.mpf(row['denominator'])
    errors[s,h]=(-1)**((s-1)//2)*(mp.zeta(s)-root)
for (s,h),err in errors.items():
    m=(s-1)//2;r=m-1;a=2*mp.pi
    c=mp.mpf(2)/factorial(2*m)*(a*a/(1+a*a))**r
    C=(1+a*a)**r
    lower=c*(mp.mpf(factorial(2*h))/factorial(2*h+r))**2
    eps=mp.mpf('1e-100')
    assert lower*errors[3,h]*(1-eps)<=err<=C*errors[3,h]*(1+eps)
print('PASS numerical corroboration:',len(rows),'archived exact-root error comparisons; no Apéry identity used')
