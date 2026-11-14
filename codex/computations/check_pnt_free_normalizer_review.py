#!/usr/bin/env python3
"""Standard-library exact checks for the independent normalizer audit."""
from fractions import Fraction as F
from math import ceil,isqrt,lcm
from pathlib import Path
import json


def primes_through(n):
    a=bytearray(b'\x01')*(n+1)
    a[:2]=b'\x00\x00'
    for p in range(2,isqrt(n)+1):
        if a[p]:a[p*p:n+1:p]=b'\x00'*((n-p*p)//p+1)
    return [p for p in range(2,n+1) if a[p]]


def factorial_vp(n,p):
    result=0
    while n:
        n//=p
        result+=n
    return result


def barnes_vp(n,p):
    result=0
    power=p
    while power<=n:
        r=n//power
        result+=r*n-power*r*(r+1)//2
        power*=p
    return result


def floor_log(n,p):
    result=0
    while n>=p:
        n//=p
        result+=1
    return result


def sum_even_floors(h,power):
    r=2*h//power
    numerator=4*h*r-power*r*(r+1)-2*(power%2)*((r+1)//2)
    assert numerator%4==0
    return numerator//4


def scalar_vp(K,p):
    h=37*K//40
    N=3*K//40
    result=2*h*factorial_vp(K,p)-12*h*factorial_vp(N,p)
    power=p
    while power<=2*h:
        result-=2*sum_even_floors(h,power)
        power*=p
    if p==2:result+=2*(h-1)
    return result


def log_interval(v,terms=32):
    v=F(v);shift=0
    while v<1:v*=2;shift-=1
    while v>=2:v/=2;shift+=1
    def base(w):
        z=(w-1)/(w+1)
        lo=2*sum((z**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
        tail=2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
        return lo,lo+tail
    lo,hi=base(v);l2,h2=base(F(2))
    return ((lo+shift*l2,hi+shift*h2) if shift>=0
            else (lo+shift*h2,hi+shift*l2))


for n in range(61):
    for p in primes_through(67):
        assert barnes_vp(n,p)==sum(factorial_vp(j,p) for j in range(n))
for h in range(1,50):
    for power in range(2,70):
        assert sum_even_floors(h,power)==sum(2*i//power for i in range(1,h))

source=Path(__file__).with_name('sparse_barnes_certificate_100.json')
data=json.loads(source.read_text())
a=list(map(F,data['coefficients']));c=list(map(F,data['frequencies']))
q=lcm(*(v.denominator for v in a))
W1=sum(abs(u)*v for u,v in zip(a,c))
C=ceil(2*W1+45)
assert q==40000 and W1==F(2053577,40000) and C==148
assert sum(u*v*v/2 for u,v in zip(a,c))==-F(111,160)

for K in [40,80,400,320000]:
    for p in primes_through(min(5*K,300)):
        L=floor_log(5*K,p)
        for cc in c:
            n=int(cc*K)
            assert abs(barnes_vp(n,p)-F(n*n,2*(p-1)))<=2*n*L
        S=scalar_vp(K,p)
        assert abs(S+F(111*K*K,160*(p-1)))<=36*K*L
        W=sum(aa*barnes_vp(int(cc*K),p) for aa,cc in zip(a,c))
        assert W-S>=-(2*W1+36)*K*L

M=40;K=200*M*M;H=K//M
assert 5*K<=H*H
for p in primes_through(H):
    assert floor_log(5*K,p)<=3*floor_log(H,p)

alpha=F(3,40);lam=F(37,40);G2=F(111,160)
upper=F(3,2)*G2+2*lam-12*alpha*lam-3*lam*lam
terms=[(aa*cc*cc/2,cc) for aa,cc in zip(a,c)]
terms += [(12*alpha*lam,alpha),(2*lam*lam,2*lam)]
for coeff,arg in terms:
    low,high=log_interval(arg)
    upper+=coeff*(high if coeff>=0 else low)
M=10**7
patch=3*C*log_interval(F(4))[1]/M
assert upper+patch<F(1379,1000)
print('PASS exact Barnes valuations, scalar floors, all-prime error bounds and patch comparison')
print('Sparse constants q=',q,'W1=',W1,'C=',C,'M=',M)
print('Certified total <1379/1000; decimal diagnostic:',float(upper+patch))
