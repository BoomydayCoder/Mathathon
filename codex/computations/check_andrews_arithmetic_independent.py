#!/usr/bin/env python3
"""Independent raising-operator / dynamic-programming / pole-normalizer audit."""
from fractions import Fraction as F
from math import comb,factorial,lcm

def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):out[i+j]+=u*v
    return out

def raise_poly(c,m,s):
    sh=[F(0)]*(2*len(c)-1)
    for k,v in enumerate(c):
        for j in range(2*k+1):sh[j]+=v*comb(2*k,j)*F(1,2)**(2*k-j)
    N=2*m+1;a=F(s,2)
    return mul(mul([a,F(1)],[F(comb(N,j))*(a+1)**(N-j) for j in range(N+1)]),sh)[1::2]

def multisum(m,n,j):
    lim=min(n,j-1)
    z=[F(comb(j+n-k,n-k)*comb(j-1,k),comb(n,k)) for k in range(lim+1)]
    state=[F(1)]+[F(0)]*lim
    for _ in range(m):
        state=[z[k]**2*sum((state[a]*comb(n+1,k-a) for a in range(k+1)),F(0)) for k in range(lim+1)]
    return sum((v*comb(n,k)*F(j,j-k) for k,v in enumerate(state)),F(0))

comparisons=norms=0
for m in range(1,5):
    for n in range(1,8):
        c=[F(1)]
        for s in reversed(range(n)):c=raise_poly(c,m,s)
        degree=m*n;lead=c[-1]
        p=[(-1)**(degree+k)*v/lead for k,v in enumerate(c)]
        C=comb((m+1)*n+m,n);dn=lcm(*range(1,n+1))
        for j in range(1,2*degree+3):
            value=sum(v*j**(2*k) for k,v in enumerate(c))/c[0]
            assert value==multisum(m,n,j)
            assert (dn**(2*m-2)*value).denominator==1
            assert value>=comb(j+n,n)**(2*m)
            if j>n:assert value<=F(j,j-n)*C*comb(j+n,n)**(2*m)
            comparisons+=1
        K=2*degree+m+1;dK=lcm(*range(1,K+1));power=2*m+1
        M=F(m*C*C*factorial(2*K)*dn**(4*m-4)*dK**power,factorial(n)**(4*m))
        assert M.denominator==1
        harmonic=F(0);constant=F(0);zcoef=F(0)
        for j in range(1,K+1):
            harmonic+=F(1,j**power)
            val=sum(v*(-j*j)**k for k,v in enumerate(p))
            residue=F(2*(-1)**(j-1)*j*j,factorial(K-j)*factorial(K+j))*val*val
            rational=-j**(2*m)*harmonic-F(1,2*m)+F(1,2*j)
            assert (M*residue*rational).denominator==1
            constant+=residue*rational;zcoef+=residue*j**(2*m)
        assert zcoef==(-1)**m
        assert (M*constant).denominator==1
        norms+=1
print(f'PASS: {comparisons} exact raising-operator versus dynamic-programming multisums and evaluation bounds; {norms} complete termwise pole-normalizer checks, including the exact zeta coefficient')
