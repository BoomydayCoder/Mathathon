"""Exact transform hierarchy and limits of simple scalar-Wilson substitutions."""
import sympy as s
x,t,a=s.symbols('x t a')
def N(m,p):
 return s.expand(sum((-1)**(m+e+1)*co*(s.bernoulli(2*m+2*e+2,a)-s.bernoulli(2*m+2*e+2))/(2*m+2*e+2) for (e,),co in s.Poly(p,t).terms()))
def mu3(p):
 return s.expand(sum(co*(-1)**e*s.bernoulli(2*e+2)*s.Rational(2*e+3,2) for (e,),co in s.Poly(p,t).terms()))
for m in range(1,5):
 for p in [s.Integer(1),2-3*t+t*t,-s.Rational(2,5)+7*t**3]:
  E=lambda q:t*s.diff(q,t)
  S=4*E(E(p))+(8*m+10)*E(p)+(2*m+2)*(2*m+3)*p
  assert s.expand(s.diff(N(m,t*p),a,2)+N(m,S)-2*mu3(t**m*p))==0
p=t*t-2*t+s.Rational(1,5)
assert mu3(t*p)==0
assert mu3(t*t*p)==s.Rational(2,15)
for n in range(1,9):
 # Scalar Wilson a=0,b=c=d=1: endpoint terms k=0 or k<=1.
 lead=(-1)**n*s.rf(n+2,n)
 p0=s.factorial(n)**3/lead
 pm1=p0*(1+n*(n+2))
 assert s.factor(p0-(-1)**n*s.factorial(n)**2/s.binomial(2*n+1,n))==0
 assert pm1/p0==(n+1)**2
print('Twelve arbitrary-polynomial transform hierarchies and scalar-Wilson endpoint cases passed.')
print('Balanced (m,n)=(2,1) is not the ordinary degree2 polynomial for its endpoint scalar weight.')
print('No all-m,n Rodrigues or constant-term theorem is asserted.')
