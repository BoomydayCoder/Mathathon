#!/usr/bin/env python3
"""Exact jets and numerical corroboration of the proved uniform asymptotic."""
from fractions import Fraction as F
from math import comb, factorial
import mpmath as mp
R=3

def J(x):return [F(x)]+[F(0)]*R
def add(a,b):return [x+y for x,y in zip(a,b)]
def neg(a):return [-x for x in a]
def sub(a,b):return add(a,neg(b))
def mul(a,b):return [sum((a[j]*b[k-j] for j in range(k+1)),F(0)) for k in range(R+1)]
def inv(a):
    assert a[0]
    b=[1/a[0]]
    for k in range(1,R+1):b.append(-sum((a[j]*b[k-j] for j in range(1,k+1)),F(0))/a[0])
    return b
def div(a,b):return mul(a,inv(b))
def scale(a,c):return [x*c for x in a]
def affine(c,d):return [F(c),F(d)]+[F(0)]*(R-1)
def alpha(n):
    k=n+1
    return scale(div(mul(affine(k*k,-1),affine(k*k,-1)),affine(3*k*k,-1)),F(k,2))
def beta(n):
    return sub(sub(scale(affine(3*n*(n+1)+1,-1),F(2*n+1)),alpha(n)),alpha(n-1))
def direct_A(n):
    a=J(0);w=J(1)
    for k in range(n+1):
        if k:w=div(w,affine(1,F(-1,k*k)))
        a=add(a,scale(w,F(comb(n,k)**2*comb(n+k,k)**2)))
    return a

def mpf(x):return mp.mpf(x.numerator)/x.denominator

def main():
    A=[J(1),direct_A(1)];b=[J(0),inv(alpha(0))]
    for n in range(1,40):
        A.append(div(sub(mul(beta(n),A[n]),mul(alpha(n-1),A[n-1])),alpha(n)))
        b.append(div(sub(mul(beta(n),b[n]),mul(alpha(n-1),b[n-1])),alpha(n)))
    for n in range(13):assert A[n]==direct_A(n)
    for n in range(1,41):
        assert sub(mul(b[n],A[n-1]),mul(b[n-1],A[n]))==inv(alpha(n-1))
    q=[div(y,x) for x,y in zip(A,b)]
    assert q[1]==add(inv(affine(1,-1)),inv(affine(5,-1)))
    mp.mp.dps=220
    expected=[(-1)**r*mp.mpf(2)**(2*r+1)*mp.pi**(2*r)/factorial(2*r+2) for r in range(R+1)]
    distances=[]
    for n in [10,20,40]:
        errors=[mp.zeta(2*r+3)-mpf(q[n][r]) for r in range(R+1)]
        assert errors[0]>0
        ratios=[e/errors[0] for e in errors]
        assert all((-1)**r*errors[r]>0 for r in range(R+1))
        distances.append([abs(x-y) for x,y in zip(ratios,expected)])
        print('n',n,'coefficient error ratios',*[mp.nstr(x,14) for x in ratios[1:]])
    assert all(distances[-1][r]<distances[0][r] for r in range(1,R+1))
    print('limits',*[mp.nstr(x,14) for x in expected[1:]])
    # Direct finite-sum A versus recurrence numerator at complex parameters.
    for sig in [mp.mpf('.25'),mp.mpf('-.25'),mp.mpc(0,'.25'),mp.mpc('.15','.15')]:
        def al(n):
            k=n+1
            return k*(k*k-sig)**2/(2*(3*k*k-sig))
        aa=[mp.mpf(1),(5-sig)/(1-sig)];bb=[mp.mpf(0),1/al(0)]
        for n in range(1,40):
            be=(2*n+1)*(3*n*(n+1)+1-sig)-al(n)-al(n-1)
            aa.append((be*aa[n]-al(n-1)*aa[n-1])/al(n))
            bb.append((be*bb[n]-al(n-1)*bb[n-1])/al(n))
        rt=mp.sqrt(sig);G=mp.pi*rt/mp.sin(mp.pi*rt)
        Z=-(mp.digamma(1+rt)+mp.digamma(1-rt)-2*mp.digamma(1))/(2*sig)
        ds=[]
        for n in [10,40]:
            w=mp.mpf(1);direct=mp.mpf(0)
            for k in range(n+1):
                if k:w/=1-sig/(k*k)
                direct+=comb(n,k)**2*comb(n+k,k)**2*w
            assert abs(direct/aa[n]-1)<mp.mpf('1e-190')
            base_error=mp.zeta(3)-mpf(q[n][0])
            ratio=(Z-bb[n]/aa[n])/base_error
            ds.append(abs(ratio-G**(-2)))
            assert abs(aa[n]/mpf(A[n][0])-1)<=mp.pi/2-1
        assert ds[-1]<ds[0]
    print('PASS: 13 independent finite-sum jets; 40 exact Casoratians; 3 fixed-order sign/limit checks; 4 complex-parameter convergence checks.')

if __name__=='__main__':main()
