#!/usr/bin/env python3
"""Exact small dictionary checks; no conjectural all-h identity is asserted."""
from fractions import Fraction as Q
from math import comb,factorial
import sympy as S
from explore_rank_one_subspaces import beta_moments
# Prevost Eq53 evaluated at y=n, Wilson normalization removed.
for n in range(16):
 total=sum(S.rf(-n,k)**2*S.rf(n+1,k)**2/S.factorial(k)**4 for k in range(n+1))
 assert total==sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))
# Rational Markov samples: coefficient of zeta is zero at m+2 nodes.
for ss in (3,5,7,9):
 m=(ss-1)//2
 for j in range(m+2,m+9):
  nodes=list(range(1,m+2))+[j]
  alpha=Q(0)
  for l in nodes:
   denominator=1
   for k in nodes:
    if k!=l:denominator*=k*k-l*l
   alpha+=Q(l**(2*m),denominator)
  assert alpha==0
# Different first approximant in Ferguson's fixed-measure family.
mu1=Q(1);mu2=1-Q(1,2**3)
assert mu1**2/mu2==Q(8,7)!=Q(5,4)
# Formal negative-parameter Wilson polynomial fails for this continuous weight.
t=S.symbols('t');a=b=c=S.Integer(1);d=S.Integer(-4)
formal=(a+b)*(a+c)*(a+d)-(a+b+c+d)*(a*a+t)
assert S.expand(formal)==t-11
K,beta=beta_moments(3,2)
assert beta[1]/beta[0]==Q(1,5)
print('PASS primary-source dictionary: all-half Wilson/Apéry, rational Markov samples, distinct Ferguson start, negative-Wilson counterexample')

# Source-guided numerator lead: finite verification only.
x,y=S.symbols('x y')
A=[S.Integer(1),S.Integer(5)];B=[S.Integer(0),S.Integer(6)]
for h in range(1,9):
 if h>=2:
  nn=h-1
  for seq in (A,B):seq.append(S.cancel(((2*nn+1)*(17*nn*nn+17*nn+5)*seq[nn]-nn**3*seq[nn-1])/(nn+1)**3))
 pi=S.expand(sum(S.binomial(h,k)*S.binomial(h+k,k)*S.prod(x-j for j in range(k))*S.prod(x+j+1 for j in range(k))/S.factorial(k)**2 for k in range(h+1)))
 divided=S.Poly(S.cancel((pi-pi.subs(x,y))/(x-y)),x)
 theta=sum(divided.nth(j)*j*S.bernoulli(j-1,0) for j in range(1,divided.degree()+1))
 C=lambda a:S.harmonic(a,3)*pi.subs(x,a)+theta.subs(y,a)/2
 assert S.simplify(C(h)-B[h])==0
 assert S.simplify(C(h-1)-(B[h]+B[h-1])/6)==0
 assert S.simplify(pi.subs(x,h-1)-(A[h]+A[h-1])/6)==0
print('PASS finite-only associated-polynomial difference lead through h=8; no all-h claim')
