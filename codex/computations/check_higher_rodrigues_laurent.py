#!/usr/bin/env python3
"""Exact algebra checks; high-precision values are corroboration, not proof."""
from fractions import Fraction
from math import comb
import sympy as s
import mpmath as mp

x, V, a, Y = s.symbols('x V a Y')
p = 2*x**3-3*x-2
value_polynomial = 256*V**3-2036001*V**2-17376*V-256
assert s.resultant(p, V*(2-x)**4-16*(1+x)**4, x) == 16*value_polynomial
assert s.discriminant(value_polynomial,V) < 0
assert value_polynomial.subs(V,s.Rational(1,256)) != 0
assert s.expand((2-x)*(1+x)**2+x**3) == 3*x+2
assert s.factor((2-x)*(1+x)**2-x**3) == -p
assert s.simplify((2*(1+x)/(2-x))**4).subs(x,-s.Rational(2,3)) == s.Rational(1,256)
curvature = 6*a/(a*a+Y)+4*(1-a)/((1-a)**2+Y)-2*(a+2)/((a+2)**2+Y)
positive = 12*((3*a+1)*Y+a*(1-a)*(a+2))/((Y+a*a)*(Y+(1-a)**2)*(Y+(a+2)**2))
assert s.factor(curvature-positive) == 0
# Reflection cancels all cotangent poles exactly: K''' sin^4/pi^4 = -4-2cos(2pi z).
z = s.symbols('z')
assert s.trigsimp(s.diff(s.pi*s.cot(s.pi*z),z,3)*s.sin(s.pi*z)**4/s.pi**4+4+2*s.cos(2*s.pi*z)) == 0
# Exact all-m root simplicity derivative, degree cancellations and value products.
for m in range(1,8):
    plus=(m-x)*(1+x)**m-x**(m+1)
    minus=(m-x)*(1+x)**m+x**(m+1)
    ratio=(m-x)*(1+x)**m/x**(m+1)
    assert s.factor(s.diff(ratio,x)/ratio+m*(m+1)/(x*(1+x)*(m-x)))==0
    assert s.Poly(plus,x).degree()==m+1
    assert s.Poly(minus,x).degree()==m-1
    assert s.gcd(plus,s.diff(plus,x))==1
    if m>1: assert s.gcd(minus,s.diff(minus,x))==1
    plus_product=(-1)**(m+1)*plus.subs(x,-1)/plus.subs(x,m)
    minus_product=(-1)**(m-1)*minus.subs(x,-1)/minus.subs(x,m)
    assert (m**(m+1)*plus_product)**(2*m)==1
    assert (m**(m-1)*minus_product)**(2*m)==s.Rational(1,m)**(4*m)
# General covariance block determinant, independently from block algebra.
for m in range(1,6):
    u = s.symbols('u',positive=True)
    aa=u/(1+u); beta=u*(m-u)/(m*(1+u)**2)
    gamma=u*(m-u)/(m*(1+u)); delta=2*u*(m-u)/m
    H = (aa*s.eye(2*m)+beta*s.ones(2*m)).row_join(gamma*s.ones(2*m,1))
    H = H.col_join((gamma*s.ones(1,2*m)).row_join(s.Matrix([[delta]])))
    expected = 2*s.Rational(m+1,m)*(m-u)*(u/(1+u))**(2*m+1)
    # Rank-one block formula avoids expensive symbolic determinant expansion.
    assert s.factor(aa**(2*m-1)*(delta*(aa+2*m*beta)-2*m*gamma**2)-expected)==0
    assert H.subs(u,s.Rational(1,2)).det() == expected.subs(u,s.Rational(1,2))
# Direct small Laurent coefficient checks independent of the closed binomial sum.
for m in (1,2):
    xx=s.symbols('x0:'+str(2*m)); zz=s.symbols('zz')
    C=s.prod(1+t for t in xx); P=s.prod(xx)
    for n in range(3):
        # z extraction first, then actual polynomial coefficient extraction in all xi.
        polynomial=s.expand(sum(comb(m*n,k)**2*C**(n+k)*P**(m*n-k) for k in range(m*n+1)))
        coefficient=s.Poly(polynomial,*xx).coeff_monomial(P**(m*n))
        expected=sum(comb(m*n,k)**2*comb(n+k,k)**(2*m) for k in range(m*n+1))
        assert coefficient==expected
print('PASS exact resultant, branch, reflection, curvature, Hessian, and constant-term checks')

def form(n):
    aa=[(-1)**k*comb(2*n,k)*comb(n+k,k)**2 for k in range(2*n+1)]
    h4=[Fraction(0)]; h5=[Fraction(0)]
    for k in range(1,2*n+1):
        h4.append(h4[-1]+Fraction(1,k**4));h5.append(h5[-1]+Fraction(1,k**5))
    A=4*sum(t*t for t in aa)
    B=4*sum((aa[k]**2*h5[k] for k in range(2*n+1)),Fraction(0))
    B-=2*sum((aa[k]*aa[l]*(h4[l]-h4[k])/(l-k) for l in range(2*n+1) for k in range(l)),Fraction(0))
    return A,B

mp.mp.dps=80
real_x=mp.findroot(lambda t:2*t**3-3*t-2,1.5)
ts=mp.findroot(lambda t:2*t**3-3*t+2,mp.mpc('.74','.37'))
psi=lambda t:6*t*mp.log(t)+4*(1-t)*mp.log(1-t)-2*(t+2)*mp.log(t+2)
psi1=lambda t:6*mp.log(t)-4*mp.log(1-t)-2*mp.log(t+2)
psi2=lambda t:6/t+4/(1-t)-2/(t+2)
f=lambda t:(1-t)**2/(t**3*(t+2))
eta=16*((1-ts)/(ts+2))**4
lam=(2*(1+real_x)/(2-real_x))**4
C=8*mp.pi**2*mp.sqrt(2*mp.pi)*f(ts)/mp.sqrt(psi2(ts))
assert abs(psi1(ts)-2j*mp.pi)<mp.mpf('1e-75')
assert abs(mp.exp(mp.log(16)+psi(ts)-2j*mp.pi*ts)-eta)<mp.mpf('1e-75')
assert abs(abs(eta)**2*lam-1)<mp.mpf('1e-75')
print('positive lambda:',mp.nstr(lam,30),'complex eta:',mp.nstr(eta,30))
print('complex coefficient C:',mp.nstr(C,30),'sharp modulus:',mp.nstr(abs(eta),30))
# At non-small phases, normalized ratios tend to 1. The bounded scaled-error
# display remains meaningful near zeros of the leading cosine.
for n in (5,10,20,40,80):
    A,B=form(n)
    with mp.workdps(12*n+100):
        ts=mp.findroot(lambda t:2*t**3-3*t+2,mp.mpc('.74','.37'))
        eta=16*((1-ts)/(ts+2))**4
        C=8*mp.pi**2*mp.sqrt(2*mp.pi)*f(ts)/mp.sqrt(psi2(ts))
        L=A*mp.zeta(5)-mp.mpf(B.numerator)/B.denominator
        leading=mp.re(C*eta**n)/(3*mp.mpf(n)**mp.mpf('2.5'))
        scaled=(L-leading)*mp.mpf(n)**mp.mpf('3.5')/abs(eta)**n
        print('n=',n,'actual=',mp.nstr(L,12),'leading=',mp.nstr(leading,12),'scaled error=',mp.nstr(scaled,12))
print('Numerics only corroborate; the contour maximum and uniform-Stirling argument are in the note.')
