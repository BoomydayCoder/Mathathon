"""Exact finite checks, not a proof of the recurrence. Python standard library."""
from fractions import Fraction as Q
from math import comb

def H(n,s=1):
    return sum((Q(1,j**s) for j in range(1,n+1)),Q(0))

def pair(n):
    A=B=Q(0)
    for k in range(n+1):
        a=comb(n,k)**2*comb(n+k,k)**2
        b2=a*(2*H(k)-H(n+k)-H(n-k))
        A+=a
        B+=a*H(k,3)+b2*H(k,2)
    return A,B

values=[pair(n) for n in range(21)]
for n in range(1,20):
    for i in (0,1):
        assert (n+1)**3*values[n+1][i]-(34*n**3+51*n*n+27*n+5)*values[n][i]+n**3*values[n-1][i]==0
for n in range(21):
    for t in (Q(1,2), Q(3,2), Q(n+2)):
        r=Q(1)
        for j in range(1,n+1): r*=Q(t-j,t+j)**2
        r/=t*t
        partial=Q(0)
        for k in range(n+1):
            a=comb(n,k)**2*comb(n+k,k)**2
            b=2*a*(2*H(k)-H(n+k)-H(n-k))
            partial+=a/(t+k)**2+b/(t+k)
        assert r==partial
print('PASS: partial fractions n=0..20 at 3 rational points; A/B recurrence n=1..19')
print(values[:5])
