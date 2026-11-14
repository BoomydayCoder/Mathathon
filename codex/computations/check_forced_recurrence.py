"""Exact finite checks of a general proved boundary-forcing formula."""
import sympy as S
from math import comb

def H(n,s=1):
 return sum((S.Rational(1,k**s) for k in range(1,n+1)),S.S.Zero)
def B(n,m):
 return sum((comb(n,k)**2*comb(n+k,k)**2*(H(k,2*m+1)+(2*H(k)-H(n+k)-H(n-k))*H(k,2*m)/m) for k in range(n+1)),S.S.Zero)
for n in range(1,25):
 E=-(2*n+1)*(3+8*n*(n+1)*(H(n-1)+H(n+1)))/(n*n*(n+1)**2)
 TB=(n+1)**3*B(n+1,2)-(2*n+1)*(17*n*n+17*n+5)*B(n,2)+n**3*B(n-1,2)
 assert S.cancel(TB+E)==0
print('PASS exact zeta5 forcing from harmonic coefficients n=1..24')
