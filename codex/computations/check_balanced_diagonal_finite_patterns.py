#!/usr/bin/env python3
"""Finite-only diagonal numerator pattern and an exact failed derivative lemma."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb,lcm
from functools import lru_cache
import json
import sympy as S
base=Path(__file__).parent
scope={}
exec((base/'check_balanced_rodrigues_integrality.py').read_text().split('rows=[];conditions=0')[0],scope)
@lru_cache(None)
def bern(k):return Q(S.bernoulli(k,0))
def moment(m,k):return (-1)**k*Q(comb(2*k+2*m+1,2*m),2*k+2)*bern(2*k+2)
rows=[];counterexamples=[]
for m in range(1,5):
 for n in range(1,16):
  cf=scope['integral_rodrigues'](m,n);den=factorial(n)**(2*m)
  c=[Q((-1)**e*v,den) for e,v in enumerate(cf)]
  a=n+1;Ra=sum(v*(-a*a)**e for e,v in enumerate(c));assert Ra.denominator==1
  quotient=[Q(0)]*(len(c)-1);acc=Q(0)
  for e in range(len(c)-1,0,-1):
   acc=c[e]-a*a*acc;quotient[e-1]=acc
  assert c[0]-a*a*acc==Ra
  alpha=a**(2*m)*Ra
  b=alpha*sum(Q(1,j**(2*m+1)) for j in range(1,a+1))+Ra*Q(1,2*m)-Ra*Q(1,2*a)-sum(v*moment(m,e) for e,v in enumerate(quotient))
  dn=lcm(*range(1,n+1));candidate=factorial(2*m)*a*dn**(2*m+1)
  assert (candidate*b).denominator==1
  L=2*m*n+2*m-1;lam=1
  for p in S.primerange(1,n+1):
   power=p
   while power*p<=L:power*=p
   lam*=power
  proved=factorial(2*m)*a*dn**(2*m-2)*lam**(2*m+1)
  assert (proved*b).denominator==1
  rows.append(dict(m=m,n=n,diagonal_numerator_denominator=str(b.denominator),candidate_multiplier_clears=True,proved_multiplier_clears=True))
  if m==2 and n in [5,6]:
   third=sum(Q(v*factorial(2*e),factorial(2*e-3)*den) for e,v in enumerate(cf) if 2*e>=3)
   failed=factorial(2*m)*dn**(2*m-2)*third
   assert failed.denominator==5
   counterexamples.append(dict(m=m,n=n,j=1,third_derivative=str(third),scaled_value=str(failed)))
print('PASS',len(rows),'finite diagonal coefficient patterns; these do not prove the sharper multiplier')
print('Exact counterexamples to proposed derivative integrality:',counterexamples)
(base/'balanced_diagonal_finite_patterns.json').write_text(json.dumps(dict(finite_only_cases=rows,failed_derivative_lemma=counterexamples),indent=2)+'\n')
