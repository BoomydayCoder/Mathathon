#!/usr/bin/env python3
"""Exact finite-product/kernel checks; numerical corroboration is labeled."""
import sympy as S
import mpmath as mp
v,q,u,y,a,b=S.symbols('v q u y a b',positive=True)
for N in range(1,7):
    finite=S.Poly(S.expand(S.prod(v+S.sin(S.pi*S.Rational(2*k-1,2*N))**2 for k in range(1,N+1))),v)
    expected=S.Poly(S.expand(S.Rational(2)**(1-2*N)*(S.chebyshevt(N,1+2*v)+1)),v)
    assert all(S.simplify(finite.nth(j)-expected.nth(j))==0 for j in range(N+1))
for N in range(1,13):
    assert S.cancel(S.chebyshevt(N,(q+1/q)/2)-(q**N+q**(-N))/2)==0
print('PASS: exact finite Chebyshev factors N1..6 and hyperbolic normalization N1..12')
h=lambda c:y/(y*y+c*c)
k=lambda c:c/(y*y+c*c)
assert S.cancel(S.diff(k(a),a)+S.diff(h(a),y))==0
assert S.cancel(k(a)*k(b)-a*b/(b*b-a*a)*(1/(y*y+a*a)-1/(y*y+b*b)))==0
integral=S.simplify(a*b/(b*b-a*a)*(S.pi/(2*a)-S.pi/(2*b)))
assert S.cancel(integral-S.pi/(2*(a+b)))==0
assert S.cancel(S.diff(integral,a,b)-S.pi/(a+b)**3)==0
print('PASS: rational kernel factorization and differentiated integral constant')
series=S.tanh(u).series(u,0,22).removeO()
for e in range(1,11):
    assert series.coeff(u,2*e-1)==2**(2*e)*(2**(2*e)-1)*S.bernoulli(2*e)/S.factorial(2*e)
    assert S.bernoulli(2*e,S.Rational(1,2))==(S.Rational(2)**(1-2*e)-1)*S.bernoulli(2*e)
print('PASS: exact tanh/Bernoulli coefficient and half-shift identities e1..10')

mp.mp.dps=65
for xx in [mp.mpf('.3'),mp.mpf(1),mp.mpf(3)]:
 for N in [1,2,3,8,17]:
    value=mp.sinh(2*xx/N)/(2*N)*sum(1/(mp.sinh(xx/N)**2+mp.sin((k-mp.mpf('.5'))*mp.pi/N)**2) for k in range(1,N+1))
    assert abs(value-mp.tanh(xx))<mp.mpf('1e-60')
for aa in [mp.mpf('.3'),mp.mpf('.5'),mp.mpf('1.5'),mp.mpf('2.5')]:
    fun=lambda yy:2*mp.pi**2*yy*mp.tanh(mp.pi*yy)*mp.sech(mp.pi*yy)**2/(yy*yy+aa*aa)
    integral=mp.quad(fun,[0,mp.mpf('.2'),1,3,mp.inf])
    assert abs(integral-2*mp.zeta(3,aa+mp.mpf('.5')))<mp.mpf('1e-55')
print('PASS numerical corroboration:15 finite-product identities and4 positive resolvents')
