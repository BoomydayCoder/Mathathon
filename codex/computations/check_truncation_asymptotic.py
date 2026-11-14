"""Exact forcing identity and high-precision checks of the truncated tail asymptotic."""
import sympy as s
import mpmath as mp
from fractions import Fraction as Q
from math import comb,factorial
n,h,u,v=s.symbols('n h u v')
S=4*n*(2*n+1)**2
S1=-4*(2*n+1)*(4*n+3)
S2=-16*(2*n+1)
d=-s.Rational(3,2)/n+1/(2*n+1)
e=s.Rational(7,4)/n**2-1/(2*n+1)**2
f=-s.Rational(15,4)/n**3+2/(2*n+1)**3
R=8*h**3+12*h*u+2*v
Rm=8*(h+d)**3+12*(h+d)*(u+e)+2*(v+f)
actual=(S*(R-Rm)+3*S1*(4*h*h+2*u)+6*S2*h+3*(n+1)*(h-1/(2*n+2)))/24
expected=-(s.Rational(287,8)*n+s.Rational(335,8)+s.Rational(25,2)/n)*h+s.Rational(351,16)+s.Rational(63,2)/n+11/n**2
assert s.cancel(actual-expected)==0
print('PASS: all-n exact forcing simplification with independent h,u,v')

print('n, E_n, ratio to (pi/2)log(n)/n^3*16^-n')
for n in [5,10,20,40,80]:
    mp.mp.dps=4*n+100
    H={}
    for q in [1,4,5]:
        H[q]=[Q(0)]
        for j in range(1,2*n+1): H[q].append(H[q][-1]+Q(1,j**q))
    a=[comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1)]
    B=sum(a[k]*H[5][k]+Q(a[k],2)*(2*H[1][k]-H[1][n+k]-H[1][n-k])*H[4][k] for k in range(n+1))
    C=sum(Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2*(H[1][n-j]+H[1][n+j]-2*H[1][j-1])/2 for j in range(1,n+1))
    E=sum(a)*mp.zeta(5)-mp.mpf((B+C).numerator)/(B+C).denominator
    leading=mp.pi/2*mp.log(n)/n**3*mp.mpf(16)**(-n)
    print(n,mp.nstr(E,8),mp.nstr(E/leading,12))
