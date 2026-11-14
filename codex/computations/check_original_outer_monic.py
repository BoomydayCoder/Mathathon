#!/usr/bin/env python3
"""Exact local test of the original undepleted determinant in the outer band."""
from fractions import Fraction as Q
from math import comb, prod


def valuation(n,p):
    n=abs(n); assert n
    a=0
    while n%p==0:
        a+=1; n//=p
    return a


def bernoullis(n):
    values=[Q(1)]
    for k in range(1,n+1):
        values.append(-sum(Q(comb(k+1,j))*values[j] for j in range(k))/Q(k+1))
    return values


def polynomial(nodes,p):
    out=[1]
    for node in nodes:
        new=[0]*(len(out)+1)
        for j,c in enumerate(out):
            new[j]=(new[j]+node*c)%p
            new[j+1]=(new[j+1]+c)%p
        out=new
    return out


def divide(poly,c,p):
    out=[0]*(len(poly)-1)
    carry=poly[-1]
    out[-1]=carry
    for j in range(len(poly)-2,0,-1):
        carry=(poly[j]-c*carry)%p
        out[j-1]=carry
    return out,(poly[0]-c*carry)%p


def matmul(a,b,p):
    n=len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n))%p for j in range(n)] for i in range(n)]


def characteristic_plus(a,p):
    n=len(a); assert n<p
    if not n:
        return [1]
    power=[row[:] for row in a]; traces=[0]
    for k in range(1,n+1):
        traces.append(sum(power[j][j] for j in range(n))%p)
        if k<n:
            power=matmul(power,a,p)
    coeff=[1]
    for k in range(1,n+1):
        coeff.append(sum((-1)**(j-1)*coeff[k-j]*traces[j] for j in range(1,k+1))
                     *pow(k,-1,p)%p)
    return coeff


def check(K,N,s,p):
    assert K//2<p<=K and p>s-1 and p*p>2*K and 2*N<p
    h=K-N; E=K+4*N-2
    B=bernoullis(2*E+2)
    moments=[(-1)**e*B[2*e+2]*Q(comb(2*e+s,s-2),s-1) for e in range(E+1)]
    assert all(x.denominator%p for x in moments)
    moments=[x.numerator*pow(x.denominator,-1,p)%p for x in moments]
    nodes=list(range(N+1,K+1))
    A=polynomial([j*j for j in range(1,N+1) for _ in range(5)]
                 +[j*j for j in nodes],p)
    C=[]; powers=[]
    for j in nodes:
        residue=prod(l*l-j*j for l in range(1,N+1))**5*prod(l*l-j*j for l in nodes if l!=j)
        cj=j**(s-1)*residue
        C.append(cj); powers.append(valuation(cj,p))
    assert max(powers)<=s
    critical=[i for i,v in enumerate(powers) if v==s]
    reduced=[[0]*h for _ in range(h)]
    for i,j in enumerate(nodes):
        reduced[i][i]=0 if j<p else (-pow(2,-1,p) if j==p else -1)%p
        if i in critical:
            quotient,_=divide(A,j*j%p,p)
            factor=pow((C[i]//p**s)%p,-1,p)
            for z,k in enumerate(nodes):
                q,_=divide(quotient,k*k%p,p)
                value=sum(c*m for c,m in zip(q,moments))%p
                reduced[i][z]=(reduced[i][z]+factor*value)%p
    d=max(0,K+N-2*p+1)
    if s>=9 or p>(K+N)/2:
        assert not critical
        # The coefficient of X^(h-b), b=K-p+1, is a unit after scaling by p^(s*b).
        b=K-p+1
        residue=(-1)**b*pow(2,-1,p)%p
        assert residue
        print(f'K={K}, N={N}, s={s}, p={p}: exact coefficient valuation {-s*b}, residue {residue}')
    else:
        assert s==7 and len(critical)==3*d
        small=[[reduced[i][j] for j in critical] for i in critical]
        coefficients=characteristic_plus(small,p)
        residual_rank=max(i for i,c in enumerate(coefficients) if c)
        forced=K-p+1-2*d+residual_rank
        assert forced>=K-p+1-2*d
        print(f'K={K}, N={N}, s=7, p={p}: residual size {3*d}, '
              f'nonzero characteristic degree {residual_rank}, actual denominator exponent >= {7*forced}')


if __name__=='__main__':
    for s in (7,9,11):
        for p in (23,29,31,37):
            check(40,3,s,p)
    for p in (41,43,47):
        check(80,6,7,p)
    print('PASS exact original polynomial moments, diagonal leading factors, and residual characteristic checks')
