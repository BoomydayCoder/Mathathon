"""Exact error-increment/bracket checks, with asymptotic diagnostics."""
import sympy as s
import mpmath as mp
from math import factorial
A=[s.Integer(1),s.Integer(5)];B=[s.Integer(0),s.Integer(6)]
for n in range(1,81):
 K=(2*n+1)*(17*n*n+17*n+5)
 A.append(s.cancel((K*A[-1]-n**3*A[-2])/(n+1)**3))
 B.append(s.cancel((K*B[-1]-n**3*B[-2])/(n+1)**3))
C={h:5*A[h]-A[h-1]for h in range(1,82)}
D={h:5*B[h]-B[h-1]for h in range(1,82)}
Q={h:s.cancel(D[h]/C[h])for h in range(1,82)}
L={h:s.Rational(h*h,12*factorial(2*h))*C[h]for h in range(1,82)}
d=lambda h:12*h*h-15*h+5
for h in range(1,82):
 assert s.cancel(Q[h]-B[h]/A[h])==s.Rational(6,h**3)/(A[h]*C[h])
 if h>=2:
  assert Q[h-1]-Q[h]==s.Rational(72*d(h),h*h*(h-1)**3)/(C[h]*C[h-1])
  W=(Q[h]-Q[h-1])*L[h]*L[h-1]
  assert W==-s.Rational(d(h),2*(h-1)*factorial(2*h)*factorial(2*h-2))
for h in [1,2,4,8]:
 total=sum(s.Rational(72*d(k),k*k*(k-1)**3)/(C[k]*C[k-1])for k in range(h+1,31))
 assert total==Q[h]-Q[30]
print('Passed81 bracket identities,80 increments/Casoratians, and4 finite tails.',flush=True)
mp.mp.dps=350
lam=17+12*mp.sqrt(2);const=4*mp.pi**3/(3+2*mp.sqrt(2))
def val(q):return mp.mpf(int(s.numer(q)))/int(s.denom(q))
for h in [10,20,40,80]:
 upper=val(Q[h])-mp.zeta(3);lower=mp.zeta(3)-val(B[h]/A[h])
 assert upper>0 and lower>0
 print('h',h,'upper/lower',mp.nstr(upper/lower,14),'normalized upper constant',mp.nstr(upper*lam**(2*h)/const,14),flush=True)
print('Asymptotic target ratio',mp.nstr(3+2*mp.sqrt(2),14),'constant',mp.nstr(const,14))
