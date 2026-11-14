"""Exact independent floor integrals and Barnes valuation checks."""
from fractions import Fraction as F
from math import floor,factorial

def J(v):k=floor(v);return k*v-F(k*(k+1),2)
def ell(u,z):return floor(u-z)+floor(u+z)+1
def B2(r):return r*r-r+F(1,6)
def frac(x):return x-floor(x)
def cuts(*us):
 return sorted({F(0),F(1,2)}|{r for u in us for r in (frac(u),1-frac(u)) if 0<r<F(1,2)})
def integ(fun,*us):
 z=cuts(*us)
 return sum((b-a)*fun((a+b)/2) for a,b in zip(z,z[1:]))
count=0
for u in [F(j,13) for j in range(46)]:
 for v in [F(j,11) for j in range(29)]:
  cov=integ(lambda z:(ell(u,z)-2*u)*(ell(v,z)-2*v),u,v)
  assert cov==(B2(frac(u-v))-B2(frac(u+v)))/2
  count+=1
for alpha in (F(1,100),F(3,40),F(1,5),F(1,2),F(4,5)):
 H=1+2*alpha
 for x in [F(j,17) for j in range(1,1200)]:
  T=floor(2*H*x);xi=H*x-F(T,2);q=floor(2*x);eta=x-F(q,2)
  direct=integ(lambda z:(T-3*ell(alpha*x,z))*(T+3*ell(alpha*x,z)-ell(x,z)-5),x,alpha*x)
  direct+=xi*(2*T-q-5)+max(0,xi-eta)
  formula=-(5-alpha)*x+J(2*H*x)/2-J(2*x)/2+J(4*alpha*x)/2-3*J((1-alpha)*x)+3*J((1+alpha)*x)-9*J(2*alpha*x)
  assert direct==formula,(alpha,x,direct,formula)

def vp(n,p):
 v=0
 while n%p==0:n//=p;v+=1
 return v
for N in range(1,51):
 B=1
 for j in range(N):B*=factorial(j)
 for p in (2,3,5,7,11,17,31,53):
  pp=p;value=F(0)
  while pp<=N:value+=pp*J(F(N,pp));pp*=p
  assert value==vp(B,p)
print('PASS',count,'covariances;5995 direct allocation integrals;400 exact Barnes valuations')
