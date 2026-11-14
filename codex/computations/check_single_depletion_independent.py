#!/usr/bin/env python3
"""Independent exact pole, moment, and singular-block audit checks."""
from fractions import Fraction as Q
from math import comb
import sympy as S

def H(n,s): return sum((Q(1,k**s) for k in range(1,n+1)),Q())
def red(x,p):
    x=Q(x)
    assert x.denominator%p
    return x.numerator*pow(x.denominator,-1,p)%p

def exceptional(L,s,p):
    C=1-L**(s-1)
    vals=[Q(2*L*C)]+[L**(s-1)*H(k,s)-H(k//L,s) for k in range(1,2*L)]
    vals += [L**(s-1)*H(L,s)-Q(L+1,2*L)]
    return any(v.numerator%p==0 or v.denominator%p==0 for v in vals)

count=0
for L in (2,3,5,7):
 for s in (7,9,53):
  for p in (23,29,41,47):
   if exceptional(L,s,p): continue
   assert p>2*L
   C=1-L**(s-1)
   # All nodes across threshold and the exceptional zero-class node.
   for j in sorted(set([1,p//L,p//L+1,p-1,p,p+1,2*p-1])):
    if j<1: continue
    c=j**(s-1)*(L**(s-1)*H(L*j,s)-H(j,s))+Q(L-1,2*L*j)
    exact=Q(p**s,j**(s-1))*c/C
    if j==p:
     expected=(L**(s-1)*H(L,s)-Q(L+1,2*L))/C
    else:
     expected=(L**(s-1)*H(L*j//p,s)-H(j//p,s))/C
    assert red(exact,p)==red(expected,p)
    assert bool(red(exact,p)) == (L*j>p)
    count+=1
print('PASS:',count,'direct exact pole residues, including s>p and the zero node')

moments=0
for L in (2,3,5):
 raw=[sum((Q(j**k,L) for j in range(L)),Q()) for k in range(25)]
 kap=[Q()]
 for k in range(1,25):
  kap.append(raw[k]-sum((comb(k-1,j-1)*kap[j]*raw[k-j] for j in range(1,k)),Q()))
 for s in (7,53):
  for e in range(12):
   bern=Q(S.bernoulli(2*e+2))
   original=(-1)**e*bern*Q(comb(2*e+s,s-2),s-1)
   depleted=(1-Q(1,L**(2*e+2)))*original
   cumulant=(-1)**e*kap[2*e+2]*Q(comb(2*e+s,s-1),L**(2*e+2))
   assert depleted==cumulant
   for p in (7,11,23):
    if L%p: assert depleted.denominator%p
   moments+=1
print('PASS:',moments,'exact cumulant/Bernoulli identities and p-integrality')

# Adversarial residual blocks: all are singular, including the zero block.
z=S.symbols('z')
for residual in (S.zeros(2),S.Matrix([[1,1],[-1,-1]]),S.Matrix([[1,2],[2,4]])):
 A=S.zeros(5)
 A[0,0]=0; A[1,1]=2; A[2,2]=3
 for i in range(2):
  for j in range(3): A[3+i,j]=(i+1)*(j+2)
  for j in range(2): A[3+i,3+j]=residual[i,j]
 f=S.Poly((z*S.eye(5)+A).det(),z)
 expected=S.expand(z*(z+2)*(z+3)*(z*S.eye(2)+residual).det())
 assert S.expand(f.as_expr()-expected)==0
 # b=4, exceptional count=2; nonzero outside roots force index>=2.
 smallest=min(k[0] for k,v in f.terms() if v)
 assert 5-smallest>=2
print('PASS: full coupled matrices with singular exceptional blocks')
alpha=Q(3,40)
cost=7*(Q(9,16)-3*alpha**2/4)
assert cost==Q(25011,6400)
assert cost-Q(1789,500)==Q(10559,32000)>Q(1,4)
print('PASS: exact PNT/Jensen constants')
