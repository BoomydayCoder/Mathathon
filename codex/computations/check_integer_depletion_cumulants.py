#!/usr/bin/env python3
"""Exact checks of the general frequency-depletion theorem; stdlib only."""
from fractions import Fraction as Q
from math import comb, gcd


def primes(n):
    out=[]
    for p in range(2,n+1):
        if all(p%q for q in out if q*q<=p): out.append(p)
    return out


def vp(n,p):
    ans=0
    while n%p==0:
        ans+=1; n//=p
    return ans


def bernoulli(N):
    out=[Q(1)]
    for n in range(1,N+1):
        out.append(-sum(Q(comb(n+1,j))*out[j] for j in range(n))/Q(n+1))
    return out


def cumulants(ell,N):
    moments=[Q(sum(j**n for j in range(ell)),ell) for n in range(N+1)]
    out=[Q(0)]
    for n in range(1,N+1):
        out.append(moments[n]-sum(comb(n-1,j-1)*out[j]*moments[n-j] for j in range(1,n)))
    return out


def binom_lucas(n,k,p):
    out=1
    while n or k:
        a,b=n%p,k%p
        if b>a:return 0
        out=out*comb(a,b)%p
        n//=p;k//=p
    return out


def main():
    B=bernoulli(80)
    comparisons=0
    for ell in range(2,21):
        kap=cumulants(ell,80)
        for s in range(3,20,2):
            C=1
            for p in primes(ell):
                if ell%p==0:C*=p**(1+vp(s-1,p))
            for k in range(2,81,2):
                moment=(ell**k-1)*B[k]*Q(comb(k+s-2,s-2),s-1)
                assert kap[k]==Q(ell**k-1,k)*B[k]
                assert moment==kap[k]*comb(k+s-2,s-1)
                assert (C*moment).denominator==1
                comparisons+=1
    sharp=0
    for s in range(3,102,2):
        for p in primes(31):
            a=vp(s-1,p)
            if p==2:k=2**a
            else:
                b=a+1
                while p**b<=s-2:b+=1
                k=p**a+(p-2)*p**b
            assert k>=2 and k%2==0 and k%(p-1)==0
            assert binom_lucas(k+s-2,s-2,p)!=0
            sharp+=1
    mixed=0
    for a in range(2,13):
        for b in range(a+1,14):
            if gcd(a,b)!=1:continue
            for s in range(3,20,2):
                for k in range(2,81,2):
                    moment=(a**k-1)*(b**k-1)*B[k]*Q(comb(k+s-2,s-2),s-1)
                    assert moment.denominator==1
                    mixed+=1
    assert (2**2-1)*(3**2-1)*B[2]*Q(comb(5,3),4)==10
    print(f'PASS: {mixed} coprime-depletion integrality checks')
    print(f'PASS: {comparisons} exact moment/cumulant/denominator checks; {sharp} Lucas sharpness constructions')


if __name__=='__main__':main()
