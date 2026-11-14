#!/usr/bin/env python3
"""Exact hyperbolic/orthogonality certificates, plus gamma corroboration."""
import sympy as S
import mpmath as mp
from math import factorial
x,t=S.symbols('x t')
D=lambda p:S.expand(-(x*x-1)*S.diff(p,x))
def homogeneous(d,j):
    coeff=[S.Integer(1)]+[S.Integer(0)]*d
    for k in range(1,j+1):
        for e in range(1,d+1): coeff[e]+=k*k*coeff[e-1]
    return coeff[d]
def a(m,j):
    return S.factorial(2*j)/2*4**(m-j)*homogeneous(m-j,j)
cur=(x-1)/2
for m in range(1,11):
    cur=D(D(cur))
    expected=sum(a(m,j)*x*(x*x-1)**j for j in range(1,m+1))
    assert S.expand(cur-expected)==0
    for j in range(1,m+1):
        previous=a(m-1,j) if j<=m-1 else 0
        left=a(m-1,j-1) if j>=2 else 0
        if m>1: assert a(m,j)==4*j*j*previous+2*j*(2*j-1)*left
    assert a(m,m)==S.factorial(2*m)/2
print('PASS exact differentiated hyperbolic expansions and complete-homogeneous coefficients m=1..10')

def moment(s,e):
    return (-1)**e*S.bernoulli(2*e+2)*S.binomial(2*e+s,s-1)/(2*e+2)
def B(k): return -S.Rational(1,2) if k==1 else S.bernoulli(k)
def bernoulli_functional(poly):
    return sum(co*B(k[0]) for k,co in S.Poly(poly,x).terms())
for degree in range(1,25):
    assert bernoulli_functional((x+1)**degree-x**degree)==(1 if degree==1 else 0)
checks=0
for m in range(1,11):
    N=2*m+1
    p=sum((-1)**(m-e)*S.binomial(N,2*e)*t**e/S.Integer(N) for e in range(m+1))
    assert S.Poly(p,t).LC()==1
    assert all(v.q==1 for v in S.Poly(N*p,t).all_coeffs())
    Q=x**N*(x+1)**N
    even=x**N*sum(S.binomial(N,2*e)*x**(2*e) for e in range(m+1))
    assert S.expand(even-(Q-Q.subs(x,x-1))/2)==0
    for j in range(1,m+1):
        actual=sum(S.Poly(p,t).nth(e)*moment(2*j+1,e+m-j) for e in range(m+1))
        alternate=(-1)**j*bernoulli_functional(S.diff(even,x,2*j-1))/(N*S.factorial(2*j))
        assert actual==alternate==0
        checks+=1
print('PASS:',checks,'exact multiple-orthogonality equations and Bernoulli finite-difference certificates')
assert moment(3,1)==moment(3,2)==S.Rational(1,12)
assert moment(5,0)==S.Rational(5,12) and moment(5,1)==S.Rational(7,24)
rho0=moment(5,0)-S.pi**2*moment(3,1)/3
rho1=moment(5,1)-S.pi**2*moment(3,2)/3
assert S.simplify(rho1/rho0-(21-2*S.pi**2)/(30-2*S.pi**2))==0
print('PASS exact s5 component masses and irrational mean expression')

mp.mp.dps=60
for m in range(1,6):
 for y in [mp.mpf('.1'),mp.mpf('.7'),mp.mpf(3)]:
    u=mp.pi*y;c=1/mp.sinh(u)**2
    density=2*mp.pi**(2*m)*y**(2*m+1)/factorial(2*m)*mp.coth(u)*sum(int(a(m,j))*c**j for j in range(1,m+1))
    gamma_sum=mp.mpf(0)
    for j in range(1,m+1):
        const=int(a(m,j))*mp.pi**(2*m-2*j)/(2*mp.pi*factorial(2*m))
        gamma_sum+=const*y**(2*(m-j))*abs(mp.gamma(1+1j*y)**(2*j+1)*mp.gamma(1j*y)/mp.gamma(2j*y))**2
    assert abs(gamma_sum/density-1)<mp.mpf('1e-55')
print('PASS numerical corroboration: all gamma component constants m=1..5 at three scales')
