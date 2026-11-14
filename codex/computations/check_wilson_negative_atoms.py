#!/usr/bin/env python3
"""Exact checks for Bernoulli continuation and negative Wilson atoms."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb

@lru_cache(None)
def bern(n):
    if n==0:return F(1)
    return -sum((F(comb(n+1,k))*bern(k) for k in range(n)),F(0))/F(n+1)
def bp(n,x):return sum((F(comb(n,k))*bern(k)*x**(n-k) for k in range(n+1)),F(0))
def moment(k,tau):
    if not tau:return (-1)**k*(2*k+1)*bp(2*k,F(1,2))
    return F((-1)**k,k+1)*(bp(2*k+2,F(1,2)+tau)-bp(2*k+2,F(1,2)))/tau**2
def atoms(tau):
    K=(tau+F(1,2)).numerator//(tau+F(1,2)).denominator
    eps=tau-K
    return eps,[( -(tau-j-F(1,2))**2,2*(tau-j-F(1,2))/tau**2) for j in range(K) if tau-j-F(1,2)]
def mul(p,q):
    r=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):r[i+j]+=a*b
    return r
def D(n):
    p=[F(1)]
    for j in range(n):p=mul(p,[(F(j)+F(1,2))**2,F(1)])
    return p
def den(n,tau):
    p=F(1)
    for j in range(1,n+1):p*=j*j-tau*tau
    return p
def mu(p,tau):return sum((x*moment(k,tau) for k,x in enumerate(p)),F(0))
def rank(a):
    a=[r[:] for r in a];row=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(row,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[row],a[pivot]=a[pivot],a[row];v=a[row][col]
        a[row]=[x/v for x in a[row]]
        for i in range(row+1,len(a)):
            v=a[i][col]
            a[i]=[x-v*y for x,y in zip(a[i],a[row])]
        row+=1
        if row==len(a):break
    return row

def main():
    taus=[F(1,2),F(2,3),F(3,4),F(1),F(5,4),F(3,2),F(2),F(5,2),F(7,3),F(7,2),F(4),F(23,7)]
    cnt=0
    for tau in taus:
        eps,parts=atoms(tau)
        assert -F(1,2)<=eps<F(1,2)
        assert all(node<0 and weight>0 for node,weight in parts)
        assert eps**2/tau**2+sum((w for _,w in parts),F(0))==1
        for k in range(13):
            rhs=eps**2/tau**2*moment(k,abs(eps))+sum((w*x**k for x,w in parts),F(0))
            assert moment(k,tau)==rhs
            cnt+=1
        for j in range(6):
            for k in range(6):
                assert mu(mul(D(j),D(k)),tau)==den(j,tau)*den(k,tau)/F(j+k+1)
        if tau>F(1,2):
            assert (1/tau)**2/(1/(2*tau-1))==(2*tau-1)/tau**2
            assert -tau*(tau-1)-F(1,4)==-(tau-F(1,2))**2
    for N in range(1,8):
        eps,parts=atoms(F(N));assert eps==0 and len(parts)==N
        assert sorted(parts)==sorted([(-(F(r)+F(1,2))**2,F(2*r+1,N*N)) for r in range(N)])
        for d in [N-1,N,N+1]:
            assert rank([[moment(j+k,F(N)) for k in range(d+1)] for j in range(d+1)])==N
    # Finite negative witnesses for parameters just beyond the positivity boundary.
    for tau,L in [(F(3,5),256),(F(2,3),30),(F(3,4),12)]:
        b=1-tau;a=[F(1)]
        for k in range(1,L+1):a.append(a[-1]*(b+k-1)/k)
        ff=mul(a,a)
        norm=sum((v/F(j+1) for j,v in enumerate(ff)),F(0))
        xnorm=sum((v/(F(j+1)*F(j+2)) for j,v in enumerate(ff)),F(0))
        deriv=[(j+1)*a[j+1] for j in range(L)]
        # integral u(1-u)^2*(D_u f)^2 = 2/[(r+2)(r+3)(r+4)] coefficientwise.
        kinetic=sum((2*v/(F(j+2)*F(j+3)*F(j+4)) for j,v in enumerate(mul(deriv,deriv))),F(0))
        assert kinetic+(1-tau*tau)*xnorm-F(1,4)*norm<0
    # All Jacobi atom functions: exact Rodrigues coefficients and L2 overlaps.
    def rising(x,n):
        z=F(1)
        for k in range(n):z*=x+k
        return z
    jac_count=0
    for tau in taus+[F(11,2),F(47,8)]:
        family=[]
        for j in range(8):
            if not j<tau-F(1,2):continue
            ss=tau-j-1;alpha=2*ss+1
            from math import factorial
            pref=rising(alpha+1,j)/factorial(j)
            v=[pref*rising(F(-j),k)*rising(j+alpha+1,k)/(rising(alpha+1,k)*factorial(k)) for k in range(j+1)]
            for k in range(j+1):
                coeff=(-k*(k-1)-(2*ss+3)*k+j*(2*tau-j))*v[k]
                if k<j:coeff+=(k+1)*(k+2*ss+2)*v[k+1]
                assert coeff==0
            assert sum((c/(ss+k+1) for k,c in enumerate(v)),F(0))==1/tau
            assert sum((c/(2*ss+k+1) for k,c in enumerate(mul(v,v))),F(0))==1/(2*ss+1)
            for old_s,old_v in family:
                assert sum((c/(ss+old_s+k+1) for k,c in enumerate(mul(v,old_v))),F(0))==0
            family.append((ss,v));jac_count+=1
    print(f'PASS: {cnt} exact continued moments; 432 Gram identities; 21 integer-rank checks; 3 finite negative witnesses; {jac_count} Jacobi atom functions and their orthogonality.')

if __name__=='__main__':main()
