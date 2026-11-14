"""Independent exact checks for the shared coefficient and tangent arguments."""
import sympy as s
from math import comb,factorial
n,j=s.symbols('n j',integer=True)
assert s.expand((2*n+2-j)*(2*n+1-j)-2*(2*n+1)*(n+1-j)-j*(j-1))==0
z,theta=s.symbols('z theta',real=True)
def P(n):return s.expand(sum(comb(n,k)*comb(n+k,k)*factorial(k)*z**(n-k) for k in range(n+1)))
def CS(n):
 v=s.expand(P(n).subs(z,s.I*theta));return s.expand(s.re(v)),s.expand(s.im(v))
for n0 in range(9):
 C,S=CS(n0);Cp,Sp=CS(n0+1)
 assert s.expand(C*Sp-Cp*S-theta**(2*n0+1))==0
 if n0:
  C0,S0=CS(n0-1)
  assert s.expand(Cp-2*(2*n0+1)*C+theta**2*C0)==0
  assert s.expand(Sp-2*(2*n0+1)*S+theta**2*S0)==0
  denom=2*(2*n0-1)
  for idx in range(n0-1,0,-1):denom=2*(2*idx-1)-theta**2/denom
  assert s.cancel(theta/denom-S/C)==0
 # Exact factorial-moment version of the shifted Legendre polynomial.
 u=s.symbols('u')
 v=s.Poly(s.expand(s.legendre(n0,1+2*u)),u)
 mapped=sum(c*factorial(k[0])*z**(n0-k[0]) for k,c in v.terms())
 assert s.expand(mapped-P(n0))==0
# Centered phase identity before invoking any integral.
C,S,a=s.symbols('C S a',real=True)
left=s.expand((C-s.I*S)*(s.cos(a)+s.I*s.sin(a))-(C+s.I*S)*(s.cos(a)-s.I*s.sin(a)))
assert s.expand(left-2*s.I*(C*s.sin(a)-S*s.cos(a)))==0
print('PASS: all-parameter recurrence identity; tangent Casoratian and finite fractions;')
print('      centered phase signs; exact factorial-moment Legendre map.')
