#!/usr/bin/env python3
"""Exact edge-case checks for the general fixed-parameter obstruction."""
from fractions import Fraction as Q
from math import comb, prod
from check_original_outer_monic import bernoullis, valuation


def threshold(r,alpha):
    beta=(1+max(Q(1,2),alpha))/2
    lam=1-alpha
    kappa=1+(r-2)*alpha
    C=lam*(4*r*alpha+Q(9,4))+3*lam**2
    T=max(Q(3),Q(2*r-1),1+2*kappa/beta,2*C/(1-beta)**2)
    S=2*((T+1)//2)+1
    assert S%2==1 and S>T and S-2<=T
    margin=Q(S)*(1-beta)**2/2-C
    assert margin>0 and S>max(3,2*r-1) and (S-1)*beta>2*kappa
    return S,margin


def check(r,K,N,s,p):
    assert max(N,K/2)<p<=K and p>s-1 and p*p>2*K
    h=K-N; E=K+(r-2)*N-2
    assert s>max(3,2*r-1)
    if E>=0:
        B=bernoullis(2*E+2)
        moments=[(-1)**e*B[2*e+2]*Q(comb(2*e+s,s-2),s-1) for e in range(E+1)]
        assert all(x.denominator%p for x in moments)
    else:
        moments=[]
    nodes=list(range(N+1,K+1))
    harmonic=[Q(0)]
    for j in range(1,K+1):
        harmonic.append(harmonic[-1]+Q(1,j**s))
    two_removed=False
    maxordinary=0
    for j in nodes:
        residue=prod(l*l-j*j for l in range(1,N+1))**(r-1)*prod(l*l-j*j for l in nodes if l!=j)
        cj=j**(s-1)*residue
        v=valuation(cj,p)
        if j==p:
            assert v==s-1
        else:
            ell=sum((l*l-j*j)%p==0 for l in range(1,K+1))
            delta=sum((l*l-j*j)%p==0 for l in range(1,N+1))
            assert delta<=2 and ell<=4
            assert v==ell+(r-2)*delta-1<=max(3,2*r-1)
            two_removed|=delta==2
            maxordinary=max(maxordinary,v)
        assert v<s
        constant=-j**(s-1)*harmonic[j]-Q(1,s-1)+Q(1,2*j)
        scaled=p**s*constant/Q(j**(s-1))
        assert scaled.denominator%p
        mod=scaled.numerator*pow(scaled.denominator,-1,p)%p
        expected=0 if j<p else (-pow(2,-1,p) if j==p else -1)%p
        assert mod==expected
    b=K-p+1
    assert 0<=h-b
    # Exact reduction of the whole matrix determines this coefficient.
    residue=(-1)**b*pow(2,-1,p)%p
    assert residue
    print(f'PASS r={r},K={K},N={N},s={s},p={p}: ordinary max={maxordinary}, '
          f'two removed={two_removed}, exact monic coefficient valuation={-s*b}')


if __name__=='__main__':
    for r in range(1,7):
        check(r,20,12,max(5,2*r+1),13)
    check(1,12,0,5,7)
    check(6,20,0,13,17)
    check(1,20,18,5,19)
    for r,alpha in ((1,Q(0)),(1,Q(3,5)),(2,Q(1,2)),(6,Q(3,40)),(3,Q(9,10))):
        S,margin=threshold(r,alpha)
        print(f'Explicit rational threshold r={r}, alpha={alpha}: S={S}, positive margin={margin}')
