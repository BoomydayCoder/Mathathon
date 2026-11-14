#!/usr/bin/env python3
"""Exact consecutive-kernel, pure-zeta coefficient and normalizer checks."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb,lcm,prod
from functools import lru_cache
import json
import mpmath as mp
base=Path(__file__).parent
scope={}
exec((base/'check_balanced_rodrigues_integrality.py').read_text().split('rows=[];conditions=0')[0],scope)
construct=scope['integral_rodrigues']

def g(a,K,x):return 1/prod((a+j)**2-x*x for j in range(K))
count=0
for n in range(1,7):
 for K in range(1,6):
  a=Q(n+1)
  @lru_cache(None)
  def diff(s,x):
   if s==0:return g(a,K,x)
   return (diff(s-1,x+Q(1,2))-diff(s-1,x-Q(1,2)))/(2*x)
  for x in [Q(1,7),Q(2,7),Q(3,7)]:
   assert diff(n,x)==prod(range(K,K+n))*g(a-Q(n,2),K+n,x)
   count+=1
print('PASS',count,'exact iterated consecutive-kernel identities')

records=[];polys={}
for m in range(1,5):
 for n in range(1,7):
  cf=construct(m,n);C=comb((m+1)*n+m,n)
  p=[Q((-1)**(m*n+e)*v,C) for e,v in enumerate(cf)]
  polys[m,n]=p
  for a in sorted(set([n//2+1,n+1])):
   K=m*n+m+1;b=a+K-1
   cs=[]
   for j in range(a,b+1):
    pj=sum(v*(-j*j)**e for e,v in enumerate(p))
    cs.append(pj/prod(ell*ell-j*j for ell in range(a,b+1) if ell!=j))
   assert sum(c*j**(2*m) for c,j in zip(cs,range(a,b+1)))==(-1)**m
   q=(-1)**m*sum(c*(j**(2*m)*sum(Q(1,k**(2*m+1)) for k in range(1,j+1))+Q(1,2*m)-Q(1,2*j)) for c,j in zip(cs,range(a,b+1)))
   factden=factorial(n)**(2*m)*factorial(a-1)**2
   assert factorial(2*b)%factden==0
   M=m*C*(factorial(2*b)//factden)*lcm(*range(1,n+1))**(2*m-2)*lcm(*range(1,b+1))**(2*m+1)
   assert (M*q).denominator==1
   assert all((M*c/Q(2*m*j*j)).denominator==1 for c,j in zip(cs,range(a,b+1)))
   records.append(dict(m=m,n=n,a=a,q_numerator=str(q.numerator),q_denominator=str(q.denominator),normalizer_digits=len(str(M))))
print('PASS',len(records),'pure-zeta coefficient and complete normalizer checks')

mp.mp.dps=80

def density(m,y):
 cs={1:1}
 for _ in range(1,m):
  nxt={}
  for j,c in cs.items():
   nxt[j]=nxt.get(j,0)+4*j*j*c
   nxt[j+1]=nxt.get(j+1,0)+2*(j+1)*(2*j+1)*c
  cs=nxt
 csch=1/mp.sinh(mp.pi*y)**2
 return 2*mp.pi**(2*m)/factorial(2*m)*y**(2*m+1)*mp.coth(mp.pi*y)*sum(c*csch**j for j,c in cs.items())

quads=[]
for m,n,a,K in [(2,2,2,3),(2,2,2,7),(2,4,3,11),(3,2,2,10)]:
 r=n//2;D=m*n+m+1;p=polys[m,n]
 pe=lambda t:sum(mp.mpf(c.numerator)/c.denominator*t**e for e,c in enumerate(p))
 S=lambda t:t*mp.fprod((t+j*j)**(2*m+2) for j in range(1,r))*(t+r*r)**(2*m+1)
 left=mp.quad(lambda y:pe(y*y)/mp.fprod(y*y+(a+j)**2 for j in range(K))*density(m,y),[0,1,3,8,mp.inf])
 right=mp.rf(K,n)/mp.rf(D,n)*mp.quad(lambda y:S(y*y)/mp.fprod(y*y+(a-r+j)**2 for j in range(K+n))*density(m,y),[0,1,3,8,mp.inf])
 assert left>0 and right>0 and abs(left-right)<mp.mpf('1e-60')*max(1,abs(left),abs(right))
 if K==D:
  row=next(v for v in records if (v['m'],v['n'],v['a'])==(m,n,a))
  q=mp.mpf(row['q_numerator'])/mp.mpf(row['q_denominator'])
  assert abs((-1)**m*(mp.zeta(2*m+1)-q)-left)<mp.mpf('1e-60')
 quads.append(dict(m=m,n=n,a=a,K=K,positive_remainder=mp.nstr(right,25)))
 print('PASS positive block transport:',m,n,a,K,mp.nstr(right,12))
(base/'balanced_consecutive_pole_checks.json').write_text(json.dumps(dict(exact_kernel_cases=count,exact_arithmetic_cases=records,positive_quadratures=quads),indent=2)+'\n')
