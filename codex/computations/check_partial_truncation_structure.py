"""Exact cutoff structure and high-precision partial-tail diagnostics."""
from fractions import Fraction as Q
from math import factorial,comb,isqrt
from functools import cache
import mpmath as mp
@cache
def H(n,s=1):
    return sum((Q(1,j**s) for j in range(1,n+1)),Q(0))
def valuation2(q):
    a,b=q.numerator,q.denominator
    out=0
    while a%2==0: a//=2;out+=1
    while b%2==0: b//=2;out-=1
    return out
for n in range(1,81):
    h=[H(n-j)+H(n+j)-2*H(j-1) for j in range(1,n+1)]
    assert all(v!=0 for v in h)
    assert all(x>y for x,y in zip(h,h[1:]))
    for j,v in enumerate(h,1):
        if 2*j>n:
            assert valuation2(v)==-((n+j).bit_length()-1)
        a=Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2
        assert a* j**4*comb(n,j)**2*comb(n+j,j)**2==1
    minimum=max(j for j,v in enumerate(h,1) if v>0)
    if n>=3:
        # Use high-precision formula only for candidate localization;
        # the exact harmonic test itself has already determined the minimum.
        mp.mp.dps=40
        t0=(mp.sqrt(1+8*n*(n+1))-1)/4
        assert minimum in (int(mp.floor(t0)),int(mp.ceil(t0)))
print('PASS exact monotonicity, nonvanishing, 2-adic valuation, reciprocal identity and cutoff candidates')
for n in (20,40,80,160):
    mp.mp.dps=4*n+100
    A=B=Q(0)
    root=[]
    for j in range(n+1):
        a=comb(n,j)**2*comb(n+j,j)**2
        b=2*a*(2*H(j)-H(n+j)-H(n-j))
        A+=a;B+=a*H(j,5)+b*H(j,4)/4
        if j:
            aj=Q(factorial(j-1)**2*factorial(n-j),factorial(n+j))**2
            root.append(aj*(H(n-j)+H(n+j)-2*H(j-1))/2)
    errors=[]
    for k in range(n+1):
        rat=B+sum(root[:k],Q(0))
        errors.append(mp.mpf(A.numerator)*mp.zeta(5)-mp.mpf(rat.numerator)/rat.denominator)
    minimum=min(range(n+1),key=lambda k:errors[k])
    best=min(range(n+1),key=lambda k:abs(errors[k]))
    assert minimum==max(j for j in range(1,n+1) if H(n-j)+H(n+j)-2*H(j-1)>0)
    assert errors[minimum]<0 and errors[0]>0 and errors[-1]>0
    crosses=[k for k in range(n) if errors[k]*errors[k+1]<0]
    assert len(crosses)==2
    center=n/mp.sqrt(2)+(mp.sqrt(2)-1)/4
    delta=mp.sqrt(n/(4*mp.sqrt(2))*(mp.log(n)/2+mp.log(mp.pi**mp.mpf('1.5')/(6*2**mp.mpf('.25')))))
    print('n',n,'minimum',minimum,'best absolute',best,'crossing intervals',crosses,
          'predicted',mp.nstr(center-delta,8),mp.nstr(center+delta,8))
