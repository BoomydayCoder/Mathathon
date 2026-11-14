#!/usr/bin/env python3
"""Exact algebra/finite normality checks; labeled analytic corroboration."""
import sympy as S
import mpmath as mp
u,v,w,Y,t=S.symbols('u v w Y t')
J={}
for k in range(1,11):
    J[k]=S.expand(sum((-1)**((l-1)//2)*S.binomial(k,l)*(w-S.Rational(1,4))**(k-l)*w**((l-1)//2) for l in range(1,k+1,2)))
    plus=(Y*Y-S.Rational(1,4)+S.I*Y)**k
    minus=(Y*Y-S.Rational(1,4)-S.I*Y)**k
    assert S.expand(plus-minus-2*S.I*Y*J[k].subs(w,Y*Y))==0
    assert S.degree(J[k],w)==k-1 and S.Poly(J[k],w).LC()==k
print('PASS: exact jump polynomials and triangular leading coefficients k=1..10')
count=0
for m in range(1,6):
 for n in range(2,7):
    P=[sum((k+1)*(r+2)*v**r for r in range(m)) for k in range(n)]
    jump=S.Poly(S.expand(sum(P[k].subs(v,-v)*J[k] for k in range(1,n))),w)
    transformed=[]
    for r in range(n-1):
        expr=S.cancel((1+u)**(m-1)*jump.nth(r).subs(v,u/(1+u)))
        poly=S.Poly(expr,u)
        assert poly.degree()<=m-1
        transformed.append(poly)
    assert not transformed[-1].is_zero
    assert S.cancel(transformed[-1].as_expr()-(n-1)*(1+u)**(m-1)*P[n-1].subs(v,-u/(1+u)))==0
    count+=1
print('PASS:',count,'exact rectangular degree-preservation and nonzero-top checks')

def moment(s,e):return (-1)**e*S.bernoulli(2*e+2)*S.binomial(2*e+s,s-1)/(2*e+2)
count=0
for m in range(1,5):
 for n in range(1,4):
    degree=m*n
    matrix=S.Matrix([[moment(2*j+1,e+k+m-j) for e in range(degree)] for j in range(1,m+1) for k in range(n)])
    rhs=S.Matrix([-moment(2*j+1,degree+k+m-j) for j in range(1,m+1) for k in range(n)])
    assert matrix.det()!=0
    sol=matrix.inv()*rhs
    p=S.Poly(t**degree+sum(sol[e]*t**e for e in range(degree)),t)
    assert S.gcd(p,p.diff()).degree()==0
    assert p.count_roots(0,S.oo)==degree
    count+=1
print('PASS:',count,'exact finite balanced systems and Sturm certificates for simple positive roots')

mp.mp.dps=65
F=lambda z:mp.asinh(1/mp.sqrt(z))**2/mp.pi**2
rho=lambda x:mp.acosh(1/mp.sqrt(x))
mass=mp.quad(rho,[0,mp.mpf('.25'),mp.mpf('.75'),1])
assert abs(mass-1)<mp.mpf('1e-55')
for x in [mp.mpf('.01'),mp.mpf('.3'),mp.mpf(2),mp.mpf(20)]:
    value=mp.quad(lambda z:rho(z)/(x+z),[0,mp.mpf('.25'),mp.mpf('.75'),1])/mp.pi**2
    assert abs(value-F(x))<mp.mpf('1e-55')
for v0 in [mp.mpf('.01'),mp.mpf('.2'),mp.mpf('.7'),mp.mpf('.99')]:
    G=F(v0/(1-v0))
    assert abs(G-rho(v0)**2/mp.pi**2)<mp.mpf('1e-55')
    upper=F(-v0+mp.mpc(0,'1e-35'))
    assert abs(upper-(G-mp.mpf('.25')-1j*mp.sqrt(G)))<mp.mpf('1e-30')
print('PASS numerical corroboration: Stieltjes mass/representation and upper-cut branches')
