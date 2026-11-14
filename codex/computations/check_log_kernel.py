"""Numerical checks of proved log-kernel asymptotics; not proofs of the limits."""
from fractions import Fraction as Q
from math import comb
import mpmath as mp

def harmonic_table(n,s):
    h=[Q(0)]
    for k in range(1,n+1):h.append(h[-1]+Q(1,k**s))
    return h

def linear(n,m):
    mp.mp.dps=3*n+100
    h=harmonic_table(2*n,1)
    hv=harmonic_table(n,2*m+1)
    hu=harmonic_table(n,2*m)
    A=0; B=Q(0)
    for k in range(n+1):
        a=comb(n,k)**2*comb(n+k,k)**2
        b=2*a*(2*h[k]-h[n+k]-h[n-k])
        A+=a; B+=a*hv[k]+b*hu[k]/(2*m)
    return A*mp.zeta(2*m+1)-mp.mpf(B.numerator)/B.denominator

print('n, J5/first_zero, (J5-first_zero)*n^8/log(n), all_m leading ratios m=2,3,4')
for n in [5,10,20,40,80]:
    mp.mp.dps=3*n+100
    h=harmonic_table(n+1,1)
    first=2*(h[n-1]+h[n+1])/Q(n*n*(n+1)**2)
    first=mp.mpf(first.numerator)/first.denominator
    J=4*linear(n,2)
    ratios=[]
    for m in [2,3,4]:
        leading=(4*mp.log(n))**(2*m-3)/(2*m*mp.factorial(2*m-3)*n**4)
        ratios.append(mp.nstr(linear(n,m)/leading,10))
    print(n,mp.nstr(J/first,14),mp.nstr((J-first)*n**8/mp.log(n),12),ratios)
