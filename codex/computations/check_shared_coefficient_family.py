"""Exact symbolic endpoint and recurrence checks; no floating-point proof."""
import sympy as S
from math import comb,factorial
z,x=S.symbols('z x')
def P(n):return S.expand(sum(comb(n,k)*comb(n+k,k)*factorial(k)*z**(n-k) for k in range(n+1)))
for n in range(9):
 f=x**n*(z-x)**n/S.factorial(n)
 at0=sum((-1)**j*S.diff(f,x,j).subs(x,0) for j in range(2*n+1))
 atz=sum((-1)**j*S.diff(f,x,j).subs(x,z) for j in range(2*n+1))
 assert S.expand(at0-(-1)**n*P(n))==0
 assert S.expand(atz-(-1)**n*P(n).subs(z,-z))==0
 for j,c in S.Poly(P(n)+P(n).subs(z,-z),z).terms():
  assert j[0]%2==0 and c.is_Integer
 if n:
  assert S.expand(P(n+1)-2*(2*n+1)*P(n)-z*z*P(n-1))==0
 Q=P(n).subs(z,-z)
 assert S.expand(Q*P(n+1)-P(n+1).subs(z,-z)*P(n)-2*z*(-z*z)**n)==0
print('PASS exact endpoints, even part, recurrence and Casoratian n=0..8')
