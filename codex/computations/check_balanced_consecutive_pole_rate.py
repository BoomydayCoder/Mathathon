#!/usr/bin/env python3
"""Exact fifth-zeta saddle algebra, plus finite positive-integral rate diagnostics."""
from pathlib import Path
from math import factorial
import json
import sympy as S
import mpmath as mp
base=Path(__file__).parent
x=S.symbols('x',real=True)
assert S.expand(S.im((x+S.I)**6*(x-6*S.I)))==2*(35*x**4-42*x**2+3)
u=(21+4*S.sqrt(21))/35
assert S.simplify(35*u*u-42*u+3)==0
assert u>1 and u<3
assert S.simplify((u+1)/(u+36)-(408+28*S.sqrt(21))/9375)==0
assert S.Rational(337500,50653)>1
print('PASS exact saddle polynomial, branch interval, closed rate and normalizer comparison')
mp.mp.dps=55

def saddle(m):
 lo=mp.mpf('0');hi=mp.mpf(m+1)
 for _ in range(240):
  v=(lo+hi)/2
  if (2*m+2)*mp.atan(1/(2*v))-mp.atan((m+1)/v)>mp.pi:lo=v
  else:hi=v
 return (lo+hi)/2

def density(m,y):
 cs={1:1}
 for _ in range(1,m):
  nxt={}
  for j,c in cs.items():
   nxt[j]=nxt.get(j,0)+4*j*j*c
   nxt[j+1]=nxt.get(j+1,0)+2*(j+1)*(2*j+1)*c
  cs=nxt
 sh=1/mp.sinh(mp.pi*y)**2
 return 2*mp.pi**(2*m)/factorial(2*m)*y**(2*m+1)*mp.coth(mp.pi*y)*sum(c*sh**j for j,c in cs.items())

rows=[]
for m in [2,3,4]:
 vstar=saddle(m);phi=(m+1)*mp.log((vstar*vstar+mp.mpf('.25'))/(vstar*vstar+(m+1)**2))
 R=2*mp.mpf(2*m+1)**(2*m+1)*mp.mpf(m+1)**(m+1)/mp.mpf(m)**m
 assert mp.log(R)+phi>3*mp.pi/8-mp.mpf(5)/16
 if m==2:assert abs(vstar*vstar*4-mp.mpf(str(S.N(u,60))))<mp.mpf('1e-50')
 for n in [4,12,30]:
  r=n//2;B=(m+1)*n+m+1
  def scaled(v):
   y=n*v;t=y*y
   logs=mp.log(t)+(2*m+2)*sum(mp.log(t+j*j) for j in range(1,r))+(2*m+1)*mp.log(t+r*r)
   logd=sum(mp.log(t+j*j) for j in range(1,B+1))
   return mp.exp(logs-logd-n*phi)*density(m,y)*n
  integral=mp.quad(scaled,[0,vstar/2,vstar,3*vstar/2,2*vstar,4*vstar,mp.inf])
  assert integral>0
  rate=phi+mp.log(integral)/n
  rows.append(dict(m=m,n=n,finite_log_error_over_n=mp.nstr(rate,25),proved_limit=mp.nstr(phi,25)))
  print('positive rate diagnostic:',m,n,mp.nstr(rate,14),'limit',mp.nstr(phi,14))
(base/'balanced_consecutive_pole_rate_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
