#!/usr/bin/env python3
"""Independent exact checks of integer depletion; no symbolic packages required."""
from fractions import Fraction as Q
from math import comb, prod, factorial, gcd


def primes(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,int(p**.5)+1))]


def vp(n,p):
    assert n
    n=abs(n); a=0
    while n%p==0:
        n//=p; a+=1
    return a


def valuation(x,p):
    return vp(x.numerator,p)-vp(x.denominator,p)


def bernoullis(n):
    b=[Q(1)]
    for k in range(1,n+1):
        b.append(-sum(Q(comb(k+1,j))*b[j] for j in range(k))/Q(k+1))
    return b


def raw_cumulants(ell,n):
    # M'=M (log M)': integer-binomial recurrence, independent of Bernoulli formula.
    raw=[Q(sum(j**k for j in range(ell)),ell) for k in range(n+1)]
    cumulants=[Q(0)]
    for k in range(1,n+1):
        cumulants.append(raw[k]-sum(Q(comb(k-1,j-1))*cumulants[j]*raw[k-j]
                                   for j in range(1,k)))
    return cumulants


def moment(s,ell,k,b):
    return (-1)**((k-2)//2)*Q(ell**k-1,s-1)*b[k]*comb(k+s-2,s-2)


def no_carries(x,y,p):
    while x or y:
        if x%p+y%p>=p:
            return False
        x//=p; y//=p
    return True


def witness(s,p):
    a=vp(s-1,p)
    if p==2:
        return 2**a
    b=a+1
    while p**b<=s-2:
        b+=1
    return p**a+(p-2)*p**b


def rational_derivative_identity(s):
    # Homogeneity reduces all a>0 to a=1. Numerator of d^(s-1)
    # [y^s/(1+y^2)] should be (s-1)! Im(1+iy)^s.
    num=[0]*s+[1]
    for r in range(s-1):
        nxt=[0]*(len(num)+1)
        for j in range(1,len(num)):
            nxt[j-1]+=j*num[j]
            nxt[j+1]+=j*num[j]
        for j in range(len(num)):
            nxt[j+1]-=2*(r+1)*num[j]
        num=nxt
    expected=[0]*len(num)
    for j in range(1,s+1,2):
        expected[j]=factorial(s-1)*(-1)**((j-1)//2)*comb(s,j)
    assert num==expected


def cyclotomic_quotient(a,b):
    numerator=[0]*(a*b+2)
    for index,value in ((a*b+1,1),(a*b,-1),(1,-1),(0,1)):
        numerator[index]+=value
    denominator=[0]*(a+b+1)
    for index,value in ((a+b,1),(a,-1),(b,-1),(0,1)):
        denominator[index]+=value
    answer=[0]*(len(numerator)-len(denominator)+1)
    for i in range(len(answer)-1,-1,-1):
        coeff=numerator[i+len(denominator)-1]
        answer[i]=coeff
        for j,value in enumerate(denominator):
            numerator[i+j]-=coeff*value
    assert not any(numerator) and sum(answer)==1
    return answer


def double_depletion_checks(bernoulli):
    checked=0
    for a in range(2,13):
        for b in range(a+1,13):
            if gcd(a,b)!=1:
                continue
            coefficients=cyclotomic_quotient(a,b)
            raw=[sum(value*j**k for j,value in enumerate(coefficients))
                 for k in range(31)]
            cumulants=[0]
            for k in range(1,31):
                cumulants.append(raw[k]-sum(comb(k-1,j-1)*cumulants[j]*raw[k-j]
                                           for j in range(1,k)))
            for k in range(2,31,2):
                assert cumulants[k]==Q((a**k-1)*(b**k-1),k)*bernoulli[k]
                for s in range(3,22,2):
                    value=moment(s,a,k,bernoulli)*(b**k-1)
                    assert value.denominator==1
                    assert value==(-1)**((k-2)//2)*cumulants[k]*comb(k+s-2,s-1)
                    checked+=1
            for n in range(1,3*a*b):
                frequency=1-int(n%a==0)-int(n%b==0)+int(n%(a*b)==0)
                assert frequency==int(n%a!=0 and n%b!=0)
    assert cyclotomic_quotient(2,3)==[1,-1,1]
    return checked


def main():
    for s in range(3,32,2):
        rational_derivative_identity(s)
    b=bernoullis(200)
    double_count=double_depletion_checks(b)
    number_moments=number_witnesses=number_actual=0
    for ell in range(2,31):
        cumulants=raw_cumulants(ell,62)
        for k in range(2,63,2):
            assert cumulants[k]==Q(ell**k-1,k)*b[k]
        ps=[p for p in primes(ell) if ell%p==0]
        for s in range(3,42,2):
            C=prod(p**(1+vp(s-1,p)) for p in ps)
            for k in range(2,63,2):
                value=moment(s,ell,k,b)
                assert value==(-1)**((k-2)//2)*cumulants[k]*comb(k+s-2,s-1)
                assert (C*value).denominator==1
                number_moments+=1
    for s in range(3,162,2):
        for p in primes(43):
            k=witness(s,p)
            assert k>=2 and k%2==0 and k%(p-1)==0
            assert no_carries(k,s-2,p)
            # Lucas assertion checked by modular binomial where sizes are modest.
            if k<=200:
                assert comb(k+s-2,s-2)%p
                for ell in (p,p*p,2*p):
                    value=moment(s,ell,k,b)
                    assert valuation(value,p)==-1-vp(s-1,p)
                    scaled=value/Q(ell**k)
                    assert valuation(scaled,p)==-k*vp(ell,p)-1-vp(s-1,p)
                    number_actual+=1
            number_witnesses+=1
            # Arbitrarily high no-carry subsequence, including p=2.
            a=vp(s-1,p)
            bb=a+1
            while p**bb<=s-2:
                bb+=1
            for bb in range(bb,bb+3):
                kk=p**a+(1 if p==2 else p-2)*p**bb
                assert kk%2==0 and kk%(p-1)==0 and no_carries(kk,s-2,p)
    print(f'PASS {number_moments} moments, {number_witnesses} Lucas witnesses, '
          f'{number_actual} actual sharp rational valuations, rescaling subsequences, '
          f'15 exact Hermite derivative identities, and {double_count} cyclotomic double-depletion moments')

if __name__=='__main__':
    main()
