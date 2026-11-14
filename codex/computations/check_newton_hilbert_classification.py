"""Exact checks for Newton/Hilbert selfadjoint classification and cubic spectra."""
import sympy as s
j,k,A,B,c,a0,l0,d,x,t,eta=s.symbols('j k A B c a0 l0 d x t eta')
a=lambda r:a0+c*r*(r+B)
lam=lambda r:l0+c*r*(r+A)
assert s.expand((a(j)-a(k))*(j+k+A)-(lam(j)-lam(k))*(j+k+B))==0
# Independently solve k=0 then k=1 for d=a_j-a_0.
res=(d-c*(B+1))*(j+1+A)-(d*(j+A)/(j+B)-c*(A+1))*(j+1+B)
assert s.factor(res-(B-A)*(d-c*j*(j+B))/(j+B))==0
print('Classification equations and two-anchor derivation: passed.',flush=True)
f=s.Function('f')(x);p=x*x*(1-x)
J=lambda q:-(x*(1-x)*s.diff(q,x,2)+(A-B*x)*s.diff(q,x))
T=lambda q:-c*(p*s.diff(q,x,2)+((A+1)*x-(B+1)*x*x)*s.diff(q,x))+(a0*x-l0)*q
anticomm=c*(x*J(f)+J(x*f))/2+(a0-c*B/2)*x*f+(c*A/2-l0)*f
assert s.expand(T(f)-anticomm)==0
h=A/(2*x)
squarepotential=s.cancel(p*h*h-(s.diff(p*h,x)+(p*h)*((A-1)/x-(B-A-1)/(1-x))))
assert s.cancel(squarepotential-(-A*A/4+A*(2*B-A)*x/4))==0
paper=2*(p*s.diff(f,x,2)+((A+1)*x-(B+1)*x*x)*s.diff(f,x)-eta*(B-eta)*x*f)
assert s.expand(T(f).subs(a0,c*eta*(B-eta))-(-c*paper/2-l0*f))==0
print('Weighted Sturm operator, anticommutator, square completion, exact primary-source map: passed.',flush=True)
for n in range(9):
 L=s.legendre(n,1-2*x)
 nn=lambda r:a0+c*r*(r+2)
 rhs=-(n+1)*nn(n)*s.legendre(n+1,1-2*x)/(2*(2*n+1))+(c*n*(n+1)/2+a0/2-l0)*L
 if n:rhs-=n*nn(n-1)*s.legendre(n-1,1-2*x)/(2*(2*n+1))
 assert s.expand(T(L).subs({A:1,B:2})-rhs)==0
print('General uniform Legendre tridiagonal action n<=8: passed.',flush=True)
Q=[s.Integer(1)]
for r in range(8):Q.append(s.expand(Q[-1]*(t+s.Rational(2*r+1,2)**2)/(r+1)**2))
mu=lambda poly:sum(co*(2*e+1)*(-1)**e*s.bernoulli(2*e,s.Rational(1,2)) for (e,),co in s.Poly(poly,t).terms())
for r in range(7):
 for rr in range(7):assert s.factor(mu(Q[r]*Q[rr])-s.Rational(1,r+rr+1))==0
for n in range(8):
 L=s.Poly(s.legendre(n,1-2*x),x)
 R=sum(co*Q[e] for (e,),co in L.terms())
 for m in range(8):
  C=lambda nn,r:s.binomial(nn,r)*s.binomial(nn+r,r)
  K=sum(C(n,r)*C(m,r) for r in range(min(n,m)+1))
  assert s.factor(R.subs(t,-s.Rational(2*m+1,2)**2)-K)==0
print('Half-grid Hilbert isometry and mixed Apéry evaluations n,m<=7: passed.',flush=True)
Q2=[s.Integer(1)]
for r in range(7):Q2.append(s.expand(Q2[-1]*(t+(r+1)**2)/(r+2)**2))
mu_tmu3=lambda poly:sum(co*(-1)**(e+1)*s.bernoulli(2*e+4)*s.Rational(2*e+5,2) for (e,),co in s.Poly(poly,t).terms())
for r in range(6):
 for rr in range(6):assert s.factor(mu_tmu3(Q2[r]*Q2[rr])-s.Rational(1,2*(r+rr+2)*(r+rr+3)))==0
print('Beta(2,2)/Wilson-one companion Gram r,k<=5: passed.',flush=True)
print('All exact checks passed. General proof and support qualifications are in the note.',flush=True)
