"""Exact root counts, confluence and certified signs; no floating-point sign tests."""
from fractions import Fraction as F
from math import comb, factorial
import sympy as S
x,t=S.symbols('x t')
for m in range(1,5):
 for n in range(1,5):
  N=m*n
  aa=[(-1)**k*comb(N,k)*comb(n+k,k)**m for k in range(N+1)]
  P=S.Poly(sum(a*x**k for k,a in enumerate(aa)),x)
  assert P.count_roots(0,1)==N
  assert S.gcd(P,P.diff()).degree()==0
  assert P.eval(1)==(-1)**N*factorial(N)//factorial(n)**m
  assert S.Rational(-aa[-2],aa[-1])==N*S.Rational(m,m+1)**m
  f=(1-x)**N
  for r in range(m):
   f=S.expand(S.diff(x**n*f,x,n)/factorial(n))
   assert S.Poly(f,x).count_roots(0,1)==(r+1)*n+(1 if r+1<m else 0)
  assert S.expand(f-P.as_expr())==0
# Generic shifted parameters: direct Mellin residues, avoiding numerical integration.
for n,m in [(1,2),(2,2),(1,3)]:
 N=n*m;alphas=[S.Rational(j+1,2*m+1) for j in range(m)]
 aa=[(-1)**k*comb(N,k)*S.prod(S.rf(k+a+1,n)/S.rf(a+1,n) for a in alphas) for k in range(N+1)]
 Mellin=sum(a/(t+k) for k,a in enumerate(aa))
 expected=(-1)**N*factorial(N)/S.prod(S.rf(a+1,n) for a in alphas)*S.prod(t-a-j for a in alphas for j in range(1,n+1))/S.prod(t+k for k in range(N+1))
 assert S.cancel(Mellin-expected)==0
# Rational Euler-transform enclosure: inner finite differences lie between zero and one.
M=80
eta=sum((sum(((-1)**j*comb(k,j)*F(1,(j+1)**5) for j in range(k+1)),F(0))/2**(k+1) for k in range(M)),F(0))
zlo=F(16,15)*eta;zhi=zlo+F(16,15*2**M)
def form(n):
 N=2*n;aa=[(-1)**k*comb(N,k)*comb(n+k,k)**2 for k in range(N+1)]
 H=[F(0)];G=[F(0)]
 for j in range(1,N+1):H.append(H[-1]+F(1,j**5));G.append(G[-1]+F(1,j**4))
 A=4*sum(a*a for a in aa)
 B=4*sum((aa[k]**2*H[k] for k in range(N+1)),F(0))-2*sum((aa[k]*aa[l]*(G[l]-G[k])/(l-k) for k in range(N+1) for l in range(k+1,N+1)),F(0))
 return A*zlo-B,A*zhi-B
lo1,hi1=form(1);lo3,hi3=form(3)
assert 0<lo1<hi1
assert lo3<hi3<0
print('PASS: 16 exact simple-root counts and endpoint/root-sum identities; shifted Mellin identities.')
print('PASS: rational Euler enclosures certify L_(1,2)>0 and L_(3,2)<0.')
print('Illustrative decimal midpoints only:',float((lo1+hi1)/2),float((lo3+hi3)/2))
