"""Exact polynomial and integer-kernel checks."""
import sympy as S
from math import comb
n=S.symbols('n');z=n*(n+1)
def phi(k):return S.prod(z-j*(j+1) for j in range(k))/S.factorial(k)**2
def L(f):return S.cancel(((n+1)**3*(f.subs(n,n+1)-f)+n**3*(f.subs(n,n-1)-f))/(2*(2*n+1)))
for k in range(12):
 assert S.cancel(L(phi(k))-k*(k+1)*phi(k)-(k*k*phi(k-1) if k else 0))==0
def K(a,b):return sum(comb(a,k)*comb(a+k,k)*comb(b,k)*comb(b+k,k) for k in range(min(a,b)+1))
for m in range(15):
 for j in range(1,25):
  assert (j+1)**3*K(j+1,m)==(2*j+1)*(j*(j+1)+2*m*(m+1)+1)*K(j,m)-j**3*K(j-1,m)
print('PASS Newton operator k=0..11 and 360 exact kernel recurrence cases')
