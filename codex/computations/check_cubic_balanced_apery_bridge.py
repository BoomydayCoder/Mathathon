#!/usr/bin/env python3
"""Independent exact original-moment and operator checks of the cubic bridge."""
from fractions import Fraction as F
from math import comb, factorial
import sympy as s

x=s.symbols('x')
for jj in range(13):
    newton_at_minus_x2=s.prod(k*k-x*x for k in range(jj))/s.factorial(jj)**2
    assert s.expand(s.rf(x,jj)*s.rf(-x,jj)-s.factorial(jj)**2*newton_at_minus_x2)==0
a,ss,j=s.symbols('a ss j')
c=ss/2
coef=(2*ss+j+3)*(c+j)**2+(c+1)**2*(2*ss+3*j+1)
assert s.expand(coef-(ss+j+1)**3)==0
# The interior coefficient induction in equation(5), after division by f_(s,n,j).
n=s.symbols('n',integer=True,positive=True)
# Ratios follow directly by cancellation of four Pochhammer endpoint products.
# f_(s+1,n-1,j) / f_(s,n,j)
r0=(n-j)*(n+2*ss+2+j)/(n*(n+2*ss+2)*(ss+j+1)**3)
# f_(s+1,n-1,j-1) / f_(s,n,j)
r1=-j/(n*(n+2*ss+2))
assert s.cancel((ss+j+1)**3*r0-(2*ss+j+2)*r1-1)==0


def leg(n):return [F((-1)**j*comb(n,j)*comb(n+j,j)) for j in range(n+1)]
def hor(c,t):return sum((v*t**k for k,v in enumerate(c)),F(0))
def harmonic(n,p):return sum((F(1,k**p) for k in range(1,n+1)),F(0))
def moment(k):return F((-1)**k)*F(s.bernoulli(2*k+2))*F(2*k+3,2)
def pairkernel(c,d):
    z=F(0);q=F(0)
    for i,ci in enumerate(c):
        for j,cj in enumerate(d):
            if i==j:z+=2*ci*cj;q-=2*ci*cj*harmonic(i,3)
            else:q+=ci*cj*(harmonic(j,2)-harmonic(i,2))/F(j-i)
    return z,q

A=[1,5];B=[F(0),F(6)]
for n in range(1,13):
    cc=(2*n+1)*(17*n*n+17*n+5)
    A.append((cc*A[-1]-n**3*A[-2])//(n+1)**3)
    B.append((cc*B[-1]-n**3*B[-2])/F((n+1)**3))

for n in range(13):
    sn=[F((-1)**j*comb(n,j)*comb(n+j+1,j)) for j in range(n+1)]
    norm=[F(0)]*(n+1);basis=[F(1)]
    for j,cc in enumerate(sn):
        assert sum((v*moment(k) for k,v in enumerate(basis)),F(0))==F(1,2*(j+1)*(j+2))
        for k,v in enumerate(basis):norm[k]+=cc*v
        nxt=[F(0)]*(len(basis)+1)
        for k,v in enumerate(basis):nxt[k]+=j*j*v/F((j+1)**2);nxt[k+1]+=v/F((j+1)**2)
        basis=nxt
    for k in range(n):
        assert sum((v*moment(k+r) for r,v in enumerate(norm)),F(0))==0
    assert sum((u*v*moment(i+j) for i,u in enumerate(norm) for j,v in enumerate(norm)),F(0))==F(n+1,4)
    aa=n+1;ra=hor(norm,-aa*aa)
    quo=[F(0)]*n
    if n:
        quo[-1]=norm[-1]
        for k in range(n-2,-1,-1):quo[k]=norm[k+1]-aa*aa*quo[k+1]
    mu=sum((v*moment(k) for k,v in enumerate(quo)),F(0))
    alpha=aa*aa*ra
    const=ra*(-aa*aa*harmonic(aa,3)-F(1,2)+F(1,2*aa))+mu
    assert alpha==F(aa*aa,6)*(A[n+1]+A[n])
    assert const==-F(aa*aa,6)*(B[n+1]+B[n])
    mixed=pairkernel(leg(n),leg(n+1))
    assert (alpha,const)==tuple(F(aa*aa,2)*v for v in mixed)

    # Directly test the exact general transform (12), not just its diagonal.
    for pole in range(1,7):
        ra=hor(norm,-pole*pole)
        quo=[F(0)]*n
        if n:
            quo[-1]=norm[-1]
            for k in range(n-2,-1,-1):quo[k]=norm[k+1]-pole*pole*quo[k+1]
        mu=sum((v*moment(k) for k,v in enumerate(quo)),F(0))
        left=(pole*pole*ra,ra*(-pole*pole*harmonic(pole,3)-F(1,2)+F(1,2*pole))+mu)
        pp=leg(pole);prev=leg(pole-1)
        pp=[v+(prev[k] if k<len(prev) else 0) for k,v in enumerate(pp)]
        z,q=pairkernel(sn,pp)
        rhs=(F(pole*pole,4)*z,F(pole*pole,4)*q-sum((v/F(2*(k+1)) for k,v in enumerate(sn)),F(0)))
        assert left==rhs,(n,pole,left,rhs)

# Rodrigues expression is independently compared as a polynomial for small s,n.
def bs(v,kk):
    cc=s.Rational(kk,2)
    AA=(x+cc)*(x+cc+1)**3
    return s.cancel((AA*v.subs(x,x+s.Rational(1,2))-AA.subs(x,-x)*v.subs(x,x-s.Rational(1,2)))/(2*x)).expand()
for ss0 in range(4):
    for nn in range(7):
        rr=s.Integer(1)
        for kk in range(ss0+nn-1,ss0-1,-1):rr=bs(rr,kk)
        predicted=s.Integer(0)
        for jj in range(nn+1):
            predicted+=s.rf(ss0+1,nn)**3*s.rf(-nn,jj)*s.rf(nn+2*ss0+2,jj)*s.rf(s.Rational(ss0,2)+x,jj)*s.rf(s.Rational(ss0,2)-x,jj)/(s.rf(ss0+1,jj)**3*s.factorial(jj))
        assert s.expand(rr-predicted)==0
# Real-model moment comparison and polynomial weighted adjoint.
tp=[F(1)]
for jj in range(31):
    ph=sum((v/F(2*(k+1)*(k+2)) for k,v in enumerate(tp)),F(0))
    assert ph==moment(jj)
    nxt=[F(0)]*(len(tp)+1)
    for k,v in enumerate(tp):nxt[k]-=k*k*v;nxt[k+1]+=(k+1)**2*v
    tp=nxt
for jj in range(13):
    f=x**jj
    adj=-x*x*(1-x)*s.diff(f,x,2)+(5*x*x-3*x)*s.diff(f,x)+(4*x-1)*f
    assert s.expand(adj-((jj+2)**2*x**(jj+1)-(jj+1)**2*x**jj))==0
# Full Newton Gram identity from original moments, independently of the model.
newtons=[];qq=[F(1)]
for jj in range(13):
    newtons.append(qq)
    nxt=[F(0)]*(len(qq)+1)
    for k,v in enumerate(qq):nxt[k]+=jj*jj*v/F((jj+1)**2);nxt[k+1]+=v/F((jj+1)**2)
    qq=nxt
for ii in range(13):
    for jj in range(13):
        val=sum((u*v*moment(k+r) for k,u in enumerate(newtons[ii]) for r,v in enumerate(newtons[jj])),F(0))
        ss0=ii+jj
        target=F(1,4) if ss0==0 else F(2*ii*jj+ss0,2*ss0*(ss0+1)*(ss0+2))
        assert val==target
print('PASS:169 exact original-measure Newton Gram entries;')
print('     31 real-operator moments and13 weighted-adjoint identities;')
print('     free-symbol raising/induction coefficients;28 Rodrigues polynomials;')
print('     13 original-moment diagonal full Apéry identities,78 general transforms,')
print('     13 Wilson norms, original orthogonality and Newton moments checked independently.')
