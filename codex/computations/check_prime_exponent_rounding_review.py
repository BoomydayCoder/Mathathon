#!/usr/bin/env python3
"""Exact checks of sign-sensitive rounding and rational polynomial content."""
from fractions import Fraction as F
from math import ceil,floor,gcd,lcm
from functools import reduce

for q in range(1,31):
    for v in range(-200,201):
        b=F(v,q)
        low=floor(b);high=ceil(b)
        assert 0<=v-q*low<=q-1
        assert 0<=q*high-v<=q-1
        for slack in range(4):
            z=ceil(-b)+slack
            assert z+b>=0 and z+low>=0


def multiply(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def power(a,q):
    r=[F(1)]
    for _ in range(q):r=multiply(r,a)
    return r


def content(a):
    denominator=lcm(*(x.denominator for x in a))
    numerator=reduce(gcd,(abs(int(x*denominator)) for x in a))
    return F(numerator,denominator)


for P in [[F(2),F(3)], [F(-3),F(0),F(2)], [F(5),F(-7),F(11)]]:
    assert content(P)==1
    for g in [F(1,6),F(17,30),F(7,2)]:
        pol=[g*x for x in P]
        assert content(pol)==g
        for q in range(1,9):
            pp=power(P,q)
            assert content(pp)==1
            for integer in [1,2,17]:
                c=F(integer)/g**q
                normalized=[c*x for x in power(pol,q)]
                assert all(x.denominator==1 for x in normalized)
                assert content(normalized)==integer
                assert c*g**q>=1
print('PASS positive and negative valuation rounding, quotient exponent bounds, and primitive-part examples')
