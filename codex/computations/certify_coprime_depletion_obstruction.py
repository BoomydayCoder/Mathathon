#!/usr/bin/env python3
"""Rational-only certificate for the uniform-trial monic determinant lower bound."""
from fractions import Fraction as F
from math import floor,ceil


def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(c,a):return (c*a[0],c*a[1]) if c>=0 else (c*a[1],c*a[0])
def point(x):return F(x),F(x)


def log12(x):
    assert 1<=x<=2
    z=(x-1)/(x+1);v=z;total=F(0)
    for k in range(40):
        total+=2*v/(2*k+1);v*=z*z
    return total,total+2*v/(81*(1-z*z))


LOG2=log12(F(2))
def log(x):
    k=0
    while x<1:x*=2;k-=1
    while x>2:x/=2;k+=1
    return add(log12(x),scale(k,LOG2))


def atan_small(x):
    assert 0<=x<=F(1,2)
    total=F(0)
    for k in range(40):total+=(-1)**k*x**(2*k+1)/(2*k+1)
    return total,total+x**81/81


PI=add(scale(16,atan_small(F(1,5))),scale(-4,atan_small(F(1,239))))
def atan(x):
    assert x>=0
    if x>1:return add(scale(F(1,2),PI),scale(-1,atan(1/x)))
    if x>F(1,2):
        return add(scale(F(1,4),PI),scale(-1,atan_small((1-x)/(1+x))))
    return atan_small(x)


def J(t,b):
    v=add(scale(t,log(b*b+t*t)),point(-3*t))
    return add(add(v,scale(t*t/b,atan(b/t))),scale(b,atan(t/b)))


def H(x):return add(scale(x*x/2,log(x)),point(-x*x/4))


a=F(3,40);lam=1-a;b=F(1,2)
V=add(scale(lam,add(add(scale(6,J(a,b)),scale(-1,J(F(1),b))),scale(-b,PI))),point(-3*lam*lam))
L=scale(5,add(add(add(H(1+a),scale(-1,H(2*a))),scale(-1,H(1-a))),point(-2*a*(1-a))))
rate=add(V,scale(-1,L))
assert rate[0]>-F(18,5)
assert 4+rate[0]>F(2,5)>F(1,3)
# Independent algebraic reduction of the field average supplied by odd-zeta agent.
field=add(add(scale(F(9,20),log(F(409,1600))),scale(-1,log(F(5,4)))),point(F(33,20)))
field=add(field,scale(F(1173,400),atan(F(3,20))))
field=add(field,scale(-F(3,2),atan(F(1,2))))
field=add(field,scale(-F(573,800),PI))
other=add(add(scale(lam,field),point(-3*lam*lam)),scale(-1,L))
assert max(rate[0],other[0])<=min(rate[1],other[1])
assert other[0]>-F(18,5)

def outward(v):
    d=10**12
    return str(F(floor(v[0]*d),d)),str(F(ceil(v[1]*d),d))

print('PASS rational-only uniform trial lower bound',outward(rate))
print('PASS V',outward(V),'L',outward(L))
print('PASS rate > -18/5; minimal s=3 arithmetic cost4 gives strict margin>2/5')
