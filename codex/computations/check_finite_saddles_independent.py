"""Exact independent verifier of saved saddle centers; no floating-point arithmetic.
Uses Fraction complex pairs, Horner/Taylor arithmetic, and 60-term atan bounds.
Does not call numerical root solvers or trust the saved modulus/margin claims.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json

def add(z,w):return z[0]+w[0],z[1]+w[1]
def mul(z,w):return z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def power(z,k):
 q=(F(1),F(0))
 for _ in range(k):q=mul(q,z)
 return q

def small(x):
 assert 0<=x<=F(1,2)
 lo=sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(60)),F(0))
 return lo,lo+x**121/F(121)
a,b=small(F(1,5)),small(F(1,239))
PI=(16*a[0]-4*b[1],16*a[1]-4*b[0])
# Exact tangent calculation behind Machin's identity.
assert (F(120,119)-F(1,239))/(1+F(120,119)*F(1,239))==1

def atan(x):
 assert x>=0
 if x>1:
  a,b=atan(1/x);return PI[0]/2-b,PI[1]/2-a
 if x>F(1,2):
  a,b=small((1-x)/(1+x));return PI[0]/4-b,PI[1]/4-a
 return small(x)

def arginterval(lo,hi):return atan(lo)[0],atan(hi)[1]
def polynomial(m,eps):
 p=[m*(comb(m,k) if k<=m else 0)-(comb(m,k-1) if 1<=k<=m+1 else 0)-(eps if k==m+1 else 0) for k in range(m+2)]
 while p[-1]==0:p.pop()
 return p

def rouche(p,z,r):
 Taylor=[]
 for k in range(len(p)):
  v=(F(0),F(0))
  for i in range(k,len(p)):v=add(v,mul((F(p[i]*comb(i,k)),F(0)),power(z,i-k)))
  Taylor.append(v)
 lhs=abs(Taylor[0][0])+abs(Taylor[0][1])
 lhs+=sum((sum(map(abs,Taylor[k]))*r**k for k in range(2,len(p))),F(0))
 rhs=max(map(abs,Taylor[1]))*r
 assert lhs<rhs
 return rhs-lhs

base=Path(__file__).resolve().parent
records=json.loads((base/'rodrigues_finite_saddles_certificate.json').read_text())['records']
assert [r['m'] for r in records]==[2,3,4,5]
safe={2:('0.0112132149','0.0112132465'),3:('0.0173831918','0.0173832381'),4:('0.0585022717','0.0585024153'),5:('0.3399158979','0.3399166785')}
checked=0
for record in records:
 m=record['m'];rows=record['components_upper_half'];assert sorted(v['j'] for v in rows)==list(range(1-m,1))
 bounds={};counts={1:0,-1:0};centers=[]
 for row in rows:
  j=row['j'];a=F(row['t_real_center']);b=F(row['t_imag_center']);r=F(row['radius']);eps=(-1)**(j+m+1)
  assert eps in (-1,1) and 0<a-r<a+r<1
  p=polynomial(m,eps);assert rouche(p,(-a,-b),r)>0
  al,ah=a-r,a+r
  if j==0:
   assert b==0;bl=bh=F(0);counts[eps]+=1
  else:
   assert b>r;bl,bh=b-r,b+r;counts[eps]+=2
  centers.append((a,b,r))
  A=arginterval(bl/ah,bh/al)
  B=arginterval(bl/(1-al),bh/(1-ah))
  D=arginterval(bl/(ah+m),bh/(al+m))
  if j:
   k=-j
   lo=(2*m+2)*A[0]+2*m*B[0]-2*D[1]
   hi=(2*m+2)*A[1]+2*m*B[1]-2*D[0]
   assert (2*k-1)*PI[1]<lo<hi<(2*k+1)*PI[0]
   phlo=2*m*(B[0]+D[0])/PI[1];phhi=2*m*(B[1]+D[1])/PI[0]
   integer=phlo.numerator//phlo.denominator
   assert integer<phlo<phhi<integer+1
  low=(m*m*((1-ah)**2+bl**2)/((ah+m)**2+bh**2))**m
  high=(m*m*((1-al)**2+bh**2)/((al+m)**2+bl**2))**m
  assert 0<low<high<1
  bounds[j]=(low,high);checked+=1
 for i,z in enumerate(centers):
  for w in centers[:i]:assert max(abs(z[0]-w[0]),abs(z[1]-w[1]))>z[2]+w[2]
 # All remaining roots: exactly one positive-real plus-branch root.
 assert counts[1]==len(polynomial(m,1))-2
 assert counts[-1]==len(polynomial(m,-1))-1
 dom=bounds[1-m]
 assert all(dom[0]>v[1] for j,v in bounds.items() if j!=1-m)
 assert F(safe[m][0])<dom[0]<dom[1]<F(safe[m][1])
 print('PASS exact root/branch/dominance certificate for m=',m)
print('PASS',checked,'upper-half/real component disks plus conjugates, all branch choices, safe decimal bounds, and strict dominance; no floats used.')
