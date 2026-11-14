"""Exact original-pencil checks for the Heun quotient and rank-two field."""
import sympy as s
from math import factorial
z,x,X=s.symbols('z x X')
def moments(M,emax):
 t=s.symbols('t');den=s.prod(t+j*j for j in range(1,M+1));out=[]
 for e in range(emax+1):
  pol,rem=s.div(t**e,den,t)
  ans=sum(c*(-1)**k*s.Rational(2*k+3,2)*s.bernoulli(2*k+2)for(k,),c in s.Poly(pol,t).terms())
  for j in range(1,M+1):
   wt=s.Rational(2*(-1)**(j-1)*j*j,factorial(M-j)*factorial(M+j))
   ans+=rem.subs(t,-j*j)*wt*(j*j*(X-s.harmonic(j,3))-s.Rational(1,2)+s.Rational(1,2*j))
  out.append(s.expand((-1)**e*ans))
 return out
for h in range(2,8):
 for q in range(1,min(h,4)+1):
  M=2*h-q+1;n=h-1;a=h-q;mm=moments(M,2*n)
  G=s.Matrix(h,h,lambda i,j:mm[i+j]);Gc=G.diff(X)
  lam=lambda k:k*(k+2-M);mu=lambda k:k*(k+1)**2*(k-M);kap=lam(n)+lam(n+1)
  R=lambda k:lam(k)+lam(k+1)-kap
  B=lambda k:lam(k)*(2*(k+1)**2-1)+mu(k)+mu(k+1)-kap*(k+1)**2
  L=lambda k:2*k*(k+1)*mu(k)
  assert all(R(k)==2*(k-n)*(k-a+1)for k in range(h))
  J=s.zeros(q)
  for k in range(a,n+1):
   J[k-a,k-a]=B(k)
   if k<n:J[k+1-a,k-a]=R(k)
   if k>a:J[k-1-a,k-a]=L(k)
  # Convert the Newton quotient to monomial quotient coordinates.
  U=s.Matrix(q,q,lambda i,j:s.Poly(s.prod(z-r*r for r in range(1,a+j+1)),z).coeff_monomial(z**(a+i)))
  Jm=U*J*U.inv()
  if a:
   Gamma=G[a:,a:]-G[a:,:a]*G[:a,:a].inv()*G[:a,a:]
  else:Gamma=G
  assert s.simplify(Gamma*Jm-Jm.T*Gamma)==s.zeros(q)
  C=Gamma.diff(X);assert C.det()!=0
  assert all(R(k)*L(k+1)>0 for k in range(a,n))
  if q==2:
   uu,vv,ww=[s.factor(w.subs(X,0))for w in [Gamma[0,0],Gamma[0,1],Gamma[1,1]]]
   S=sum(j*j for j in range(1,M+1));assert C==s.Matrix([[0,1],[1,S]])
   Delta=8*h*h-12*h+5;center=s.factor((uu*S-2*vv)/2);width=s.factor(h*h*uu/2)
   assert s.factor(Gamma.det()+(X-center)**2-width*width*Delta)==0
   assert s.factor((J.trace())**2-4*J.det()-4*h**4*Delta)==0
  print('h',h,'rank',q,'original pencil/quotient identities passed.',flush=True)
# First applicable rational Pell pair, independently from the original moments.
h=11;M=21;mm=moments(M,20);G=s.Matrix(h,h,lambda i,j:mm[i+j]);a=h-2
Gamma=G[a:,a:]-G[a:,:a]*G[:a,:a].inv()*G[:a,a:];u=Gamma[0,0];v=Gamma[0,1].subs(X,0)
center=s.factor((u*sum(j*j for j in range(1,M+1))-2*v)/2);width=s.factor(h*h*u/2)
lo=s.factor(center-29*width);hi=s.factor(center+29*width)
assert lo==s.Rational(13780864457900933987428453,11464402743063221545440000)
assert hi==s.Rational(8054563655047033145573,6700650887554237152000)
assert s.factor(Gamma.det()+(X-lo)*(X-hi))==0
print('h11 rational Pell pair passed exactly:',lo,hi,flush=True)
hs=[1,11]
for _ in range(3):hs.append(34*hs[-1]-hs[-2]-24)
for h in hs:assert s.sqrt(8*h*h-12*h+5).is_Integer
print('Pell indices checked:',hs)
