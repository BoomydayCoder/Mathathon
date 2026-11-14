"""Independent Fraction certificate: direct two-variable Jensen primitive and lead."""
from fractions import Fraction as F
from math import prod,floor,ceil

def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(q,a):return (q*a[0],q*a[1]) if q>=0 else (q*a[1],q*a[0])
def point(q):return q,q

def log(q):
 assert q>0
 z=(q-1)/(q+1)
 partial=2*sum((z**(2*k+1)/F(2*k+1) for k in range(200)),F(0))
 error=2*abs(z)**401/(401*(1-z*z))
 return partial-error,partial+error

def atan_small(q):
 assert 0<=q<=F(1,2)
 lo=sum(((-1)**k*q**(2*k+1)/F(2*k+1) for k in range(60)),F(0))
 return lo,lo+q**121/121
PI=add(scale(16,atan_small(F(1,5))),scale(-4,atan_small(F(1,239))))
def atan(q):
 if q>1:return add(scale(F(1,2),PI),scale(-1,atan(1/q)))
 if q>F(1,2):return add(scale(F(1,4),PI),scale(-1,atan_small((1-q)/(1+q))))
 return atan_small(q)

def rectangle_average(a,b):
 return add(add(scale(a,log(a*a+b*b)),point(-3*a)),add(scale(a*a/b,atan(b/a)),scale(b,atan(a/b))))
def double_log_primitive(z):
 return add(scale(z*z/2,log(z)),point(-3*z*z/4))
alpha=F(3,40);lam=1-alpha;b=F(1,2)
external=scale(lam,add(add(scale(6,rectangle_average(alpha,b)),scale(-1,rectangle_average(F(1),b))),scale(-b,PI)))
cross=add(add(double_log_primitive(1+alpha),scale(-1,double_log_primitive(1-alpha))),scale(-1,double_log_primitive(2*alpha)))
rate=add(add(point(-3*lam*lam),external),scale(-5,cross))
assert rate[0]>-F(18,5)
assert 4+rate[0]>F(2,5)>F(1,3)
# Exact cancellation of Vandermonde denominators, including its sign.
cases=0
for K in range(1,13):
 for N in range(K):
  nodes=list(range(N+1,K+1));h=len(nodes)
  Vand2=prod((r*r-j*j)**2 for i,r in enumerate(nodes) for j in nodes[:i])
  rho=[]
  for r in nodes:
   numerator=prod(j*j-r*r for j in range(1,N+1))**6
   denominator=prod(j*j-r*r for j in range(1,K+1) if j!=r)
   rho.append(F(numerator,denominator))
  for s in [3,5,7]:
   C=(2**(s-1)-1)*(3**(s-1)-1)
   direct=Vand2*prod(rho)*prod(C*r**(s-1) for r in nodes)
   expected=(-1)**(h*(h-1)//2)*C**h*prod(r**(s-1)*prod(j*j-r*r for j in range(1,N+1))**5 for r in nodes)
   assert direct==expected and direct.denominator==1
   # Complete scaling degree in the restricted Andreief lower bound.
   power=2*h*(h-1)+(10*N-2*h)*h+h+s*h
   assert power==10*N*h+(s-1)*h
   cases+=1
# Exact prime-density integral of (1-u/6) from2 to6.
assert (F(6)-F(36,12))-(F(2)-F(4,12))==F(4,3)
d=10**12
print('PASS independent exact Jensen interval:',F(floor(rate[0]*d),d),F(ceil(rate[1]*d),d))
print('PASS',cases,'leading-coefficient/sign/scale cases and exact4s/3 prime-density integral.')
