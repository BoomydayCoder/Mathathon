"""Independent high-precision consistency checks for the truncated asymptotic.
Exact rational arithmetic is used before evaluating zeta(5); the theorem is
proved by the recurrence squeeze in truncated-derivative.md.
"""
from fractions import Fraction as Q
from math import factorial,comb
from functools import cache
import mpmath as mp
mp.mp.dps=180
@cache
def H(k,r=1):
    return sum((Q(1,j**r) for j in range(1,k+1)),Q(0))
print('n, normalized truncated error; predicted limit pi/2 =',mp.nstr(mp.pi/2,25))
for n in (2,5,10,20,40):
    A=B=C=Q(0)
    for k in range(n+1):
        a=comb(n,k)**2*comb(n+k,k)**2
        b=2*a*(2*H(k)-H(n+k)-H(n-k))
        A+=a
        B+=a*H(k,5)+b*H(k,4)/4
    for j in range(1,n+1):
        C+=Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2*(H(n-j)+H(n+j)-2*H(j-1))/2
    rational=B+C
    E=mp.mpf(A.numerator)*mp.zeta(5)-mp.mpf(rational.numerator)/rational.denominator
    scale=mp.power(16,-n)*mp.log(n)/n**3
    print(n,mp.nstr(E/scale,25))
    assert E>0
