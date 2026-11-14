"""Independent cotangent-polynomial and exact central-saddle checks."""
import sympy as S
from math import comb,factorial
c,a,Y,m,t=S.symbols('c a Y m t')
for mm in range(1,7):
 r=2*mm-1
 row=[sum((-1)**j*comb(r+1,j)*(k+1-j)**r for j in range(k+2)) for k in range(r)]
 assert all(x>0 for x in row) and sum(row)==factorial(r)
 P=c
 for _ in range(r):P=S.expand(-(1+c*c)*S.diff(P,c))
 w=(c+S.I)/(c-S.I)
 rhs=sum(v*w**(k+1-mm) for k,v in enumerate(row))
 assert S.cancel(-P/(1+c*c)**mm-rhs)==0
curvature=(2*m+2)*a/(a*a+Y)+2*m*(1-a)/((1-a)**2+Y)-2*(a+m)/((a+m)**2+Y)
claimed=2*m*(m+1)*((3*a+m-1)*Y+a*(1-a)*(a+m))/((Y+a*a)*(Y+(1-a)**2)*(Y+(a+m)**2))
assert S.factor(curvature-claimed)==0
for mm,t0,eta0,C2 in [(1,1/S.sqrt(2),17-12*S.sqrt(2),2**S.Rational(3,2)*S.pi**3/(1+S.sqrt(2))**4),(2,S.Rational(2,3),S.Rational(1,256),S.pi**5/8)]:
 f=(1-t0)**mm/(t0**(mm+1)*(t0+mm))
 second=(2*mm+2)/t0+2*mm/(1-t0)-2/(t0+mm)
 eta=(mm*(1-t0)/(t0+mm))**(2*mm)
 assert S.simplify(eta-eta0)==0
 assert S.simplify(mm**2*(2*S.pi)**(2*mm+1)*f*f/second-C2)==0
# Check the beta-factor identity by exact cancellation before invoking reflection.
z=S.symbols('z')
for nn in [1,2,3]:
 for mm in [1,2,3,4]:
  kappa=S.factorial(mm*nn)/S.factorial(nn)**mm
  H=S.gamma(z)**(2*mm+2)*S.gamma(nn+1-z)**(2*mm)/S.gamma(z+mm*nn+1)**2
  beta1=S.gamma(z)*S.gamma(nn+1-z)/S.factorial(nn)
  beta2=S.gamma(z)*S.factorial(mm*nn)/S.gamma(z+mm*nn+1)
  assert S.cancel(kappa*kappa*H-beta1**(2*mm)*beta2**2)==0
print('PASS: six independent cotangent/Eulerian identities, general curvature, exact Apéry and m2 central constants, and positive beta-factor identity.')
