"""Independent exact original-moment tests for the rank-two Apéry formula."""
import sympy as s
from math import factorial
h=s.symbols('h');v=s.symbols('v');Delta=8*h*h-12*h+5;d=12*h*h-15*h+5
y=lambda k:k*k*(k+1)/2
old=(2*h-1)**2
assert s.factor((1-old)*y(h-1)+old*y(h)-h*h*d/2)==0
ys=h*h*(2*h*h-2*h+1-v)/(2*(2*h-1))
alpha=h*h*(6*h*h-8*h+3-(2*h-1)*v)/2
beta=-h*h*(3*h-2+v)/(2*(2*h-1))
assert s.factor(alpha-old*ys-(1-old)*y(h-1))==0
assert s.factor(beta-ys+y(h))==0
a0=h*(6*h*h-8*h+3);a1=h*(2*h-1)
b0=-2*(h-1)**2*(3*h-2);b1=2*(h-1)**2
assert s.factor(a0*b1-a1*b0-2*h*(h-1)**2*d)==0
A=[s.Integer(1),s.Integer(5)];B=[s.Integer(0),s.Integer(6)]
for n in range(1,13):
 K=(2*n+1)*(17*n*n+17*n+5)
 A.append(s.cancel((K*A[-1]-n**3*A[-2])/(n+1)**3))
 B.append(s.cancel((K*B[-1]-n**3*B[-2])/(n+1)**3))
X=s.symbols('X')
for hh in range(2,13):
 M=2*hh-1;mom=[]
 for e in range(2*hh-1):
  mom.append(sum(s.Rational(2*(-1)**(j-1)*j*j,factorial(M-j)*factorial(M+j))*j**(2*e)*(j*j*(X-s.harmonic(j,3))-s.Rational(1,2)+s.Rational(1,2*j))for j in range(1,M+1)))
 G=s.Matrix(hh,hh,lambda i,j:mom[i+j]);aa=hh-2
 Gam=G if not aa else G[aa:,aa:]-G[aa:,:aa]*G[:aa,:aa].inv()*G[:aa,aa:]
 u=Gam[0,0];vv=Gam[0,1].subs(X,0);S=sum(j*j for j in range(1,M+1))
 center=s.factor((u*S-2*vv)/2);width=s.factor(hh*hh*u/2);dd=8*hh*hh-12*hh+5
 C=5*A[hh]-A[hh-1];Cp=5*A[hh-1]-A[hh-2]
 D=5*B[hh]-B[hh-1];Dp=5*B[hh-1]-B[hh-2]
 for sig in [-1,1]:
  av=hh*(6*hh*hh-8*hh+3-sig*(2*hh-1)*s.sqrt(dd))
  bv=-2*(hh-1)**2*(3*hh-2+sig*s.sqrt(dd))
  got=s.radsimp((av*D+bv*Dp)/(av*C+bv*Cp))
  assert s.simplify(got-center-sig*width*s.sqrt(dd))==0
 U=hh*(6*hh*hh-8*hh+3)*C-2*(hh-1)**2*(3*hh-2)*Cp
 V=hh*(2*hh-1)*C+2*(hh-1)**2*Cp
 assert u==-s.Rational(288*(12*hh*hh-15*hh+5)**2,hh**3*(hh-1))/(U*U-dd*V*V)
 print('h',hh,'both original roots and residual norm match Apéry formulas.',flush=True)
print('All symbolic endpoint/bilinear and11 original-pencil formula tests passed.')
