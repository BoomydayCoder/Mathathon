"""Exact membership, independence, and endpoint-pairing checks for all ranks."""
import sympy as s
from math import factorial
from explore_rank_one_subspaces import beta_moments
z,x,X=s.symbols('z x X');R={};Q={};L={}
for j in range(1,10):
 _,mom=beta_moments(3,j)
 G=s.Matrix(j,j,lambda r,c:s.Rational(mom[r+c].numerator,mom[r+c].denominator))
 vv=[s.Integer(1)]if j==1 else list(-G[:-1,:-1].inv()*G[:-1,-1])+[s.Integer(1)]
 p=s.expand((-1)**(j-1)*sum(c*(-z)**r for r,c in enumerate(vv)))
 R[j]=s.expand(p/p.subs(z,0));L[j]=s.Poly(R[j],z).LC()
 Q[j]=G[-1,-1]if j==1 else s.cancel(G[-1,-1]-(G[-1,:-1]*G[:-1,:-1].inv()*G[:-1,-1])[0])
count=0
for h in range(2,10):
 for q in range(2,min(h,6)+1):
  M=2*h-q+1;a=h-q;l=a+1;n=h-1;eps=(-1)**(q-1);sgn=(-1)**q
  F={j:R[j]if 2*j>=M else s.expand(s.prod(r*r-z for r in range(2*j+1,M+1))*R[j])for j in range(l,h+1)}
  B=s.Matrix(h,q,lambda r,c:s.Poly(F[l+c],z).coeff_monomial(z**r))
  assert B.rank()==q
  def ell(p):
   assert s.degree(p,z)<=M-1
   return s.factor(sum(s.Rational(2*(-1)**(j-1)*j*j,factorial(M-j)*factorial(M+j))*p.subs(z,j*j)*(j*j*(X-s.harmonic(j,3))-s.Rational(1,2)+s.Rational(1,2*j))for j in range(1,M+1)))
  for j in range(l,h+1):
   for k in range(a):assert ell(F[j]*z**k)==0
   expected=sgn*L[h]*(X-Q[h])if j==h else sgn*eps*L[l]*(X-Q[l])if j==l else 0
   assert s.factor(ell(F[j]*z**a)-expected)==0
   count+=1
  # Check that the explicit Heun operator preserves the proposed complement.
  lam=lambda k:k*(k+2-M);kap=lam(n)+lam(n+1)
  av=(x+1)**3*(x-M)/(2*x*(2*x+1))
  def D(p):return s.cancel(av*(p.subs(x,x+1)-p)+av.subs(x,-x)*(p.subs(x,x-1)-p))
  HC=[]
  for j in range(l,h+1):
   p=F[j].subs(z,x*x);hp=s.cancel(x*x*D(p)+D(x*x*p)-D(p)-kap*x*x*p)
   hp=s.Poly(hp,x);assert all(k%2==0 for(k,),c in hp.terms())
   hc=s.Matrix([hp.coeff_monomial(x**(2*r))for r in range(h)])
   coeff=B.gauss_jordan_solve(hc)[0]
   assert B*coeff==hc
   HC.append(coeff)
  HB=s.Matrix.hstack(*HC)
  assert HB.shape==(q,q)
  print('h',h,'rank',q,'basis, all lower orthogonalities, extreme pairings, and Heun invariance passed.',flush=True)
print('Passed',count,'basis-vector endpoint identities across all tested ranks.')
