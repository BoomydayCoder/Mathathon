"""Exact checks of both interpolation operators and the zeta2 diagonal closure."""
import sympy as S
from math import comb
n=S.symbols('n');z=n*(n+1)
def phi(k):return S.prod(z-j*(j+1) for j in range(k))/S.factorial(k)**2
def psi(k):return S.prod(n-j for j in range(k))/S.factorial(k)
def B(f):return S.cancel(((n+1)**2*f.subs(n,n+1)-n*n*f.subs(n,n-1))/(2*n+1))
def Hdiff(f):return S.expand((n+1)**2*(f.subs(n,n+1)-f)+n*n*(f.subs(n,n-1)-f))
for k in range(10):
 assert S.cancel(B(phi(k))-(2*k+1)*phi(k)-(2*k*phi(k-1) if k else 0))==0
 assert S.cancel(Hdiff(psi(k))-k*(k+1)*psi(k)-(k*k*psi(k-1) if k else 0))==0
def H(a,b):return sum(comb(a,k)*comb(a+k,k)*comb(b,k) for k in range(min(a,b)+1))
for N in range(1,41):
 weights=[comb(N,k)**2*comb(N+k-1,k) for k in range(N+1)]
 assert sum((k*k+N*k-N*N)*v for k,v in enumerate(weights))==0
 C=H(N,N);Cm=H(N-1,N-1);U=H(N,N-1);V=H(N-1,N)
 assert C==2*U+V and Cm==2*V-U
 assert (N+1)**2*H(N+1,N+1)==(11*N*N+11*N+3)*C+N*N*Cm
for a in range(1,15):
 for b in range(15):
  assert (a+1)**2*H(a+1,b)==(2*a+1)*(2*b+1)*H(a,b)+a*a*H(a-1,b)
x,N=S.symbols('x N',nonzero=True)
for r in range(1,9):
 for s in range(1,9):
  degree=S.Poly((N-x)**r*(N+x)**s-x**(r+s),x).degree()
  assert (degree<=2)==((r,s) in [(1,1),(2,1),(2,2)])
print('PASS exact basis, 40 balance/closure cases, mixed recurrence, and degree classification')
