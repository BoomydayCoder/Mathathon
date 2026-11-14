#!/usr/bin/env python3
"""Symbolic certificates plus independent exact rational parts of mixed sums."""
import sympy as S
x,y,t,n,m=S.symbols('x y t n m')
N=n*(n+1);M=m*(m+1);dn=(t-n)*(t+n+1)
left=-(n+1)**3*(t-n-1)*(t-n)-(2*n+1)*(N+2*M+1)*dn-n**3*(t+n)*(t+n+1)
right=-2*(2*n+1)*(t**4-dn*(t-m-1)*(t+m))
assert S.expand(left-right)==0
assert S.expand((t-n)**2+(t+n)**2+6*(t-n)*(t+n)-4*(2*t*t-n*n))==0
assert S.expand(4*(t**4-(t*t-n*n)**2)/n**2-4*(2*t*t-n*n))==0
assert S.expand(6*(2*n+1)*(3*n*(n+1)+1)-(n+1)**3-n**3-(2*n+1)*(17*n*n+17*n+5))==0
print('PASS symbolic mixed recurrence, adjacent telescoping certificates and diagonal 17/5 reduction')

def bernoulli_functional(poly):
 p=S.Poly(poly,x)
 return S.expand(sum(p.nth(j)*S.bernoulli(j,0) for j in range(p.degree()+1))) if p else S.Integer(0)
def L(poly):return bernoulli_functional(S.diff(poly,x))
def phi(k,z):return S.prod((z-j)*(z+j+1) for j in range(k))/S.factorial(k)**2
def Pi(k,z):return S.expand(sum(phi(j,k)*phi(j,z) for j in range(k+1)))
def P(k):return S.prod(t-j for j in range(1,k+1))/S.prod(t+j for j in range(k+1))
def theta(k):return L(S.cancel((Pi(k,x)-Pi(k,y))/(x-y)))
def C(k,a):return S.harmonic(a,3)*Pi(k,a)+theta(k).subs(y,a)/2
for k in range(9):
 f=x*x*S.prod(x*x-j*j for j in range(1,k+1))/S.factorial(k+1)**2
 assert S.expand(f.subs(x,x+1)-f-(2*x+1)*phi(k,x)/(k+1))==0
 assert L(phi(k,x))==0
 assert L(x*phi(k,x))==S.Rational((-1)**k,k+1)
 assert L(x*Pi(k,x))==(1 if k==0 else 0)
 if k:
  assert S.simplify(6*C(k,k-1)-C(k,k)-C(k-1,k-1))==0
print('PASS exact Bernoulli telescopes, associated functionals and mixed balance through n=8')

# Read the rational part independently from finite pole jets of P_n P_m.
for nn in range(6):
 for mm in range(6):
  mixed=(-1)**(nn+mm)*P(nn)*P(mm)
  aa=bb=cc=S.Integer(0)
  for k in range(max(nn,mm)+1):
   jet=S.cancel((t+k)**2*mixed)
   a=jet.subs(t,-k);b=S.diff(jet,t).subs(t,-k)
   aa+=a;bb+=b;cc+=a*S.harmonic(k,3)+b*S.harmonic(k,2)/2
  assert S.simplify(aa-Pi(nn,mm))==0
  assert S.simplify(bb)==0
  assert S.simplify(cc-C(nn,mm))==0
  cert=S.cancel(mixed*(t-mm-1)*(t+mm))
  assert S.diff(cert,t).subs(t,1)==(1 if nn==0 else 0)
print('PASS independent pole-jet rational parts and exact summation boundaries for 36 mixed pairs')
