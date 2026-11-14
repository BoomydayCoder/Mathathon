#!/usr/bin/env python3
"""Full p-adic characteristic-polynomial checks of the wider-band theorem.

Uses only exact Python arithmetic, including Newton identities with explicit
precision loss when the characteristic index is divisible by p.
"""
from fractions import Fraction as Q
from math import comb, prod
from check_original_outer_monic import bernoullis, valuation, polynomial, divide


def rational_mod(x,p,modulus):
    assert x.denominator%p
    return x.numerator*pow(x.denominator,-1,modulus)%modulus


def multiply(a,b,modulus):
    n=len(a)
    columns=list(zip(*b))
    return [[sum(x*y for x,y in zip(row,column))%modulus for column in columns] for row in a]


def characteristic_mod(a,p,precision):
    n=len(a)
    padding=sum(valuation(j,p) for j in range(1,n+1))
    full=precision+padding
    modulus=p**full
    power=[row[:] for row in a]
    traces=[0]
    for k in range(1,n+1):
        traces.append(sum(power[i][i] for i in range(n))%modulus)
        if k<n:
            power=multiply(power,a,modulus)
    coefficients=[1]
    available=full
    for k in range(1,n+1):
        modulus=p**available
        numerator=sum((-1)**(j-1)*coefficients[k-j]*traces[j] for j in range(1,k+1))%modulus
        loss=valuation(k,p)
        assert numerator%p**loss==0
        available-=loss
        coefficients.append((numerator//p**loss)*pow(k//p**loss,-1,p**available)%p**available)
    assert available==precision
    return [c%p**precision for c in coefficients]


def multiplicity(coefficients,root,p):
    poly=list(reversed(coefficients))
    while len(poly)>1 and poly[-1]==0:
        poly.pop()
    count=0
    while len(poly)>1:
        quotient,remainder=divide(poly,(-root)%p,p)
        if remainder:
            break
        poly=quotient; count+=1
    return count


def check(K,N,s,p):
    h=K-N; q=K//p; d=max(0,2*q-3); E0=q*d
    assert K/4<p<=K and s>=13 and p>s-1 and p*p>2*K and 2*N<p
    precision=d*h+E0+2
    padding=sum(valuation(j,p) for j in range(1,h+1))
    modulus=p**(precision+padding)
    maxdegree=K+4*N-2
    bernoulli=bernoullis(2*maxdegree+2)
    moments=[rational_mod((-1)**e*bernoulli[2*e+2]*Q(comb(2*e+s,s-2),s-1),p,modulus)
             for e in range(maxdegree+1)]
    nodes=list(range(N+1,K+1))
    numerator=polynomial([j*j for j in range(1,N+1) for _ in range(5)]
                         +[j*j for j in nodes],modulus)
    harmonics=[Q(0)]
    for j in range(1,K+1):
        harmonics.append(harmonics[-1]+Q(1,j**s))
    scaled=[]
    for j in nodes:
        residue=prod(l*l-j*j for l in range(1,N+1))**5*prod(l*l-j*j for l in nodes if l!=j)
        cj=j**(s-1)*residue
        v=valuation(cj,p)
        assert v<=s+d
        factor=p**(s+d-v)*pow((cj//p**v)%modulus,-1,modulus)%modulus
        once,_=divide(numerator,j*j%modulus,modulus)
        row=[]
        for k in nodes:
            twice,_=divide(once,k*k%modulus,modulus)
            value=factor*sum(c*m for c,m in zip(twice,moments))%modulus
            if j==k:
                constant=-j**(s-1)*harmonics[j]-Q(1,s-1)+Q(1,2*j)
                value=(value+rational_mod(p**(s+d)*constant/Q(j**(s-1)),p,modulus))%modulus
            row.append(value)
        scaled.append(row)
    # scaled = p^d A. Transform det(ZI+p^d A) back to Q=p^E0 det(ZI+A).
    char=characteristic_mod(scaled,p,precision)
    qcoeff=[]
    for r,c in enumerate(char):
        exponent=E0-d*r
        if exponent>=0:
            value=c*p**exponent
        else:
            assert c%p**(-exponent)==0
            value=c//p**(-exponent)
        qcoeff.append(value%p**(E0+1))
    v=min(valuation(c,p) for c in qcoeff if c)
    assert v<=E0
    reduced=[c//p**v%p for c in qcoeff]
    distinct=set()
    for a in range(1,q+1):
        H=sum(Q(1,j**s) for j in range(1,a+1))
        root=rational_mod(H,p,p)
        assert root and root not in distinct
        distinct.add(root)
        count=sum(1 for j in nodes if j//p==a and j%p)
        actual=multiplicity(reduced,root,p)
        assert actual>=max(0,count-v)
    largest=max(r for r,c in enumerate(reduced) if c)
    forced=s*max(0,K-p+1-q-q*E0)
    actual_denominator=E0+s*largest-v
    assert actual_denominator>=forced
    print(f'PASS K={K}, s={s}, p={p}: repair E={E0}, Gauss content={v}, '
          f'full coefficient denominator exponent >= {actual_denominator}, theorem >= {forced}')


if __name__=='__main__':
    check(40,3,13,13)
    check(40,3,13,17)
    check(80,6,13,23)
