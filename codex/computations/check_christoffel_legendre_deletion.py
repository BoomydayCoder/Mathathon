#!/usr/bin/env python3
"""Exact finite checks of the Legendre-only cubic-tail projection."""
from fractions import Fraction as F
from math import comb


def leg(n):
    return [F((-1)**j*comb(n,j)*comb(n+j,j)) for j in range(n+1)]


def pairing(a,b):
    return sum((x*y/F((i+j+1)*(i+j+2))
                for i,x in enumerate(a) for j,y in enumerate(b)), F(0))


previous_correction=[]
removed=[]
harmonic=F(0)
for n in range(1,41):
    ln,prev=leg(n),leg(n-1)
    numerator=[ln[j]+(prev[j] if j<len(prev) else 0) for j in range(n+1)]
    q=[]
    carry=F(0)
    for j in range(n):
        carry+=numerator[j]
        q.append(carry)
    assert numerator[n]==-q[-1]
    assert q[-1]==-ln[-1]
    assert pairing(q,q)==F(2,n)
    for r in range(n-1):
        mono=[F(0)]*r+[F(1)]
        assert pairing(q,mono)==0

    correction=[sum((ln[j]/j for j in range(r+1,n+1)),F(0))
                for r in range(n)]
    assert all(correction[j]-(previous_correction[j] if j<n-1 else 0)
               ==-q[j]/n for j in range(n))
    removed.append(F(0))
    for j in range(n):
        removed[j]+=q[j]/n
    assert correction==[-c for c in removed]
    harmonic+=F(1,n**3)
    assert pairing(correction,correction)==2*harmonic
    logpair=sum((c/F((j+1)**2) for j,c in enumerate(correction)),F(0))
    assert logpair==-2*harmonic
    for r in range(n):
        mono=[F(0)]*r+[F(1)]
        assert F(1,(r+1)**2)+pairing(correction,mono)==0
    previous_correction=correction
print('PASS:40 adjacent-Legendre norms, 780 lower-degree orthogonalities,')
print('     40 exact deletion identities and Pythagorean cubic corrections,')
print('     820 exact logarithmic projection orthogonalities.')
