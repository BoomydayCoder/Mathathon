"""Finite-grid and symbolic certificates for the restricted real boundary proof."""
import sympy as s
from math import factorial
x,X,M,j=s.symbols('x X M j')
A=lambda x:(x+1)**3*(x-M)/(2*x*(2*x+1))
# Factorial ratios are replaced by their elementary adjacent ratio.
ratio_alpha=-(j+1)**4*(M-j)/(j**4*(M+j+1))
assert s.factor(A(j)-ratio_alpha*A(-j-1))==0
# Divide delta_j by alpha_j to avoid symbolic factorial simplification.
delta_over_alpha=-(M-j)/(4*j**4)
assert s.factor(A(j)*(2*j+1)/(2*j**3*(j+1)**3)-delta_over_alpha)==0
H=s.symbols('H')
rj=X-H-s.Rational(1,2)/j**2+s.Rational(1,2)/j**3
rnext=X-H-1/(j+1)**3-s.Rational(1,2)/(j+1)**2+s.Rational(1,2)/(j+1)**3
assert s.factor(rj-rnext-(2*j+1)/(2*j**3*(j+1)**3))==0
count=0
for m in range(1,13):
 aa=A(x).subs(M,m)
 def D(p):return s.cancel(aa*(p.subs(x,x+1)-p)+aa.subs(x,-x)*(p.subs(x,x-1)-p))
 weights=[s.Rational(2*(-1)**(r-1)*r*r,factorial(m-r)*factorial(m+r))for r in range(1,m+1)]
 data=[r*r*(X-s.harmonic(r,3))-s.Rational(1,2)+s.Rational(1,2*r)for r in range(1,m+1)]
 for r in range(m):
  for ss in range(m-r):
   f=x**(2*r)+s.Rational(2,5);g=x**(2*ss)-s.Rational(3,7)
   K=s.expand(f*g.subs(x,x+1)-g*f.subs(x,x+1))
   skew=s.expand(f*D(g)-g*D(f))
   lhs=s.factor(sum(w*v*skew.subs(x,k)for k,(w,v)in enumerate(zip(weights,data),1)))
   rhs=-s.Rational(m,2*factorial(m)**2)*K.subs(x,0)
   assert s.factor(lhs-rhs)==0
   dd=lambda k:s.Rational((-1)**k,2*factorial(m-k-1)*factorial(m+k)) if k>=0 else s.Rational((-1)**(-k),2*factorial(m-k-1)*factorial(m+k))
   assert sum(dd(k)*K.subs(x,k)for k in range(-m,m))==0
   count+=1
 print('M',m,'all permitted degree pairs passed.',flush=True)
print('Passed',count,'finite-grid boundary checks plus symbolic balance/source identities.')
