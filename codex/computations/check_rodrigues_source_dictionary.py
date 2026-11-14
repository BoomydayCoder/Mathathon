#!/usr/bin/env python3
"""Exact normalization checks for the primary-source parameter dictionary."""
from math import factorial,comb
from fractions import Fraction as F
import sympy as s
z=s.symbols('z')
for m in range(1,5):
 for n in range(1,5):
  M=2*m+2;a=[n+1]*M;b=[1]*(2*m)+[(m+1)*n+2]*2
  gamma=F(1)
  for aa,bb in zip(a,b):
   gamma*=F(1,factorial(aa-bb)) if aa>=bb else factorial(bb-aa-1)
  kappa=F(factorial(m*n),factorial(n)**m)
  assert gamma==kappa*kappa
  assert sum(b)-sum(a)==2
  # Independent gamma-brick products evaluated at a rational nonpole.
  ss=F(3,2); tt=ss+n+1;source=gamma
  for aa,bb in zip(a,b):
   if aa>=bb:
    for j in range(bb,aa):source*=ss+j
   else:
    for j in range(aa,bb):source/=ss+j
  target=kappa*kappa
  for j in range(1,n+1):target*=(tt-j)**(2*m)
  for j in range(m*n+1):target/=(tt+j)**2
  assert source==target
  # Multiple-integral theorem gamma prefactor.
  pref=F(factorial(m*n)**2,factorial(n)**(2*m))
  assert pref==kappa*kappa
  # Confluent monic polynomial normalization, including Mellin variable shift.
  LC=(-1)**(m*n)*comb((m+1)*n,n)**m
  rising=factorial((m+1)*n)//factorial(m*n)
  CII=F((-1)**(m*n)*factorial(m*n),rising**m)
  assert LC*CII==kappa
# The finite Fourier coefficients in Nesterenko Lemma4 match the normalized
# cotangent derivative. This tests the principal part required by that lemma.
x=s.symbols('x')
for m in range(1,5):
 row=[1]
 for r in range(2,2*m):row=[(k+1)*(row[k] if k<len(row) else 0)+(r-k)*(row[k-1] if k else 0) for k in range(r)]
 # Pair the symmetric Fourier terms to keep the exact Taylor check small.
 trig=sum(s.Rational(v,factorial(2*m-1))*s.exp(2*s.I*(k+1-m)*x) for k,v in enumerate(row))
 polynomial=s.series(trig,x,0,2*m).removeO()
 sinfactor=s.series((x/s.sin(x))**(2*m),x,0,2*m).removeO()
 assert s.series(polynomial*sinfactor,x,0,2*m).removeO().expand()==1
print('PASS exact Nesterenko brick/derivative dictionary, central integral prefactor, confluent normalization and Fourier coefficients')
