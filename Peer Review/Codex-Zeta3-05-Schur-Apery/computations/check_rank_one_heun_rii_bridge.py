"""Exact certificates for the original Schur/Apéry Heun-RII proof.

The all-h proof is in rank-one-original-apery-heun-rii-proof.md.
"""
import sympy as s
from explore_rank_one_subspaces import beta_moments
x,z,h,k=s.symbols('x z h k')
n=h-1;M=2*h
kap=1-2*n*n
lam=lambda k:k*(k+2-M)
mu=lambda k:k*(k+1)**2*(k-M)
bdiag=lambda k:s.expand(lam(k)*(2*(k+1)**2-1)+mu(k)+mu(k+1)-kap*(k+1)**2)
eta=1-2*h-2*h**3*(h+1)
assert s.factor(lam(k)+lam(k+1)-kap-2*(k-n)*(k-n+1))==0
assert s.factor(bdiag(n)-eta)==0
assert s.factor((1-M-eta)/(2*M)-h*h*(h+1)/2)==0
ratio=(h-1)*(h**3+2*h*h+3*h-2)/(8*(2*h-1))
assert s.factor((eta-bdiag(0))/(2*mu(1)*2)-ratio)==0
y=lambda q:q*q*(q+1)/2
zz=lambda q:q*q*(q+1)**2*(3*q**3+3*q-2)/(16*(2*q-1))
assert s.factor(y(h)*(1+3*ratio)-zz(h))==0
assert s.factor((eta+kap+4*(1-M)/3)/(16*(1-M)/3)-zz(h)/y(h))==0
d=lambda q:12*q*q-15*q+5
PP=204*h**5+51*h**4-134*h**3-43*h*h+30*h+12
aa=PP/(h*h*(h+1)**2*(2*h+1)*d(h))
bb=(48*h**4-12*h**3-34*h*h+3*h+7)/((h+1)*(2*h+1)*d(h))
cc=(h-1)*d(h+1)/(4*h*h*(h+1)*(2*h-1)*(2*h+1)*d(h))
FF=lambda v:(v-(2*h-1)**2)*(v-4*h*h)
assert s.factor(bb-cc*FF(0)-1)==0
assert s.factor((aa+bb)*y(h)-cc*FF(1)*y(h-1)-y(h+1))==0
assert s.factor((4*aa+bb)*zz(h)-cc*FF(4)*zz(h-1)-zz(h+1))==0
mat=s.Matrix([[0,1,-FF(0)],[y(h),y(h),-FF(1)*y(h-1)],[4*zz(h),zz(h),-FF(4)*zz(h-1)]])
det=s.factor(mat.det());print('Endpoint system determinant:',det,flush=True)
assert det!=0
alpha=2*PP/((h+1)**3*d(h));gamma=(h-1)**3*d(h+1)/(h*(h+1)**2*d(h))
assert s.factor(alpha-aa*h*h*(2*h+2)*(2*h+1)/(h+1)**2)==0
assert s.factor(gamma-cc*(h-1)**2*(2*h-1)*(2*h)*(2*h+1)*(2*h+2)/(h+1)**2)==0
K=lambda r:(2*r+1)*(17*r*r+17*r+5)
u,v=s.symbols('u v')
up=(K(h)*u-h**3*v)/(h+1)**3
umm=(K(h-1)*v-h**3*u)/(h-1)**3
assert s.factor(5*up-u-alpha*(5*u-v)+gamma*(5*v-umm))==0
assert s.factor(25*(h-1)**3-5*K(h-1)+h**3+12*h*d(h))==0
print('All symbolic-in-h endpoint, recurrence, and Apéry-elimination identities passed.',flush=True)
A=[s.Integer(1),s.Integer(5)];B=[s.Integer(0),s.Integer(6)]
for j in range(1,11):
 A.append(s.cancel((K(j)*A[-1]-j**3*A[-2])/(j+1)**3))
 B.append(s.cancel((K(j)*B[-1]-j**3*B[-2])/(j+1)**3))
RR={};SS={}
for hh in range(1,11):
 _,mom=beta_moments(3,hh)
 G=s.Matrix(hh,hh,lambda i,j:s.Rational(mom[i+j].numerator,mom[i+j].denominator))
 vv=[s.Integer(1)]if hh==1 else list(-G[:-1,:-1].inv()*G[:-1,-1])+[s.Integer(1)]
 Q=s.expand((-1)**(hh-1)*sum(c*(-z)**j for j,c in enumerate(vv)))
 q=G[-1,-1] if hh==1 else s.cancel(G[-1,-1]-(G[-1,:-1]*G[:-1,:-1].inv()*G[:-1,-1])[0])
 R=s.expand(Q/Q.subs(z,0));RR[hh]=R
 beta=lambda j:-j*j*s.harmonic(j,3)-s.Rational(1,2)+s.Rational(1,2*j)
 S=s.interpolate([(j*j,beta(j)*R.subs(z,j*j))for j in range(1,hh+2)],z)
 assert s.degree(S,z)<=hh
 for j in range(1,2*hh+1):assert s.factor(S.subs(z,j*j)-beta(j)*R.subs(z,j*j))==0
 assert s.factor(s.Poly(S,z).LC()+q*s.Poly(R,z).LC())==0
 SS[hh]=S
 def Dop(p):
  a=(x+1)**3*(x-2*hh)/(2*x*(2*x+1))
  return s.cancel(a*(p.subs(x,x+1)-p)+a.subs(x,-x)*(p.subs(x,x-1)-p))
 p=R.subs(z,x*x);H=x*x*Dop(p)+Dop(x*x*p)-Dop(p)-(1-2*(hh-1)**2)*x*x*p
 assert s.cancel(H-eta.subs(h,hh)*p)==0
 assert R.subs(z,1)==y(s.Integer(hh)) and R.subs(z,4)==zz(s.Integer(hh))
 scale=s.Rational(hh*hh,12*s.factorial(2*hh))
 assert s.Poly(R,z).LC()==scale*(5*A[hh]-A[hh-1])
 assert q==(5*B[hh]-B[hh-1])/(5*A[hh]-A[hh-1])
 if hh>=3:
  j=hh-1;av=aa.subs(h,j);bv=bb.subs(h,j);cv=cc.subs(h,j);fac=(z-(2*j-1)**2)*(z-4*j*j)
  assert s.expand(R-(av*z+bv)*RR[j]+cv*fac*RR[j-1])==0
  assert s.expand(S-(av*z+bv)*SS[j]+cv*fac*SS[j-1])==0
 print('h',hh,'original Gram, Heun equation, endpoints, RII pair, and Apéry formula passed.',flush=True)
print('All exact bridge checks passed. Independent proof audit remains separate.')
