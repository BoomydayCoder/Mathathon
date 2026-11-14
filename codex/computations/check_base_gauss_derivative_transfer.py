#!/usr/bin/env python3
"""Exact base-Gauss derivative transfer and denominator checks."""
import sympy as S
import mpmath as mp
t,z,a=S.symbols('t z a',positive=True)
Theta=lambda f:S.expand(t*S.diff(f,t))
for m in range(1,6):
    value=(z+t)**(-m)
    for j in range(1,m):value=S.cancel(2*t*S.diff(value,t)+(2*j+1)*value)
    value=S.cancel(value/(m*S.factorial2(2*m-1)))
    other=2*S.diff(1/(a*a+t),a,2*m-2)/S.factorial(2*m)
    assert S.cancel(value.subs(z,a*a)-other)==0
print('PASS: exact base-kernel Euler/parameter derivative identities m1..5')

def moment(m,k):return (-1)**k*S.bernoulli(2*k,S.Rational(1,2))*S.prod(2*k+2*j+1 for j in range(m))/(m*S.factorial2(2*m-1))
for m in range(1,9):
 for n in range(1,9):
    left=S.binomial(4*n+2*m-1,2*m-2)*moment(1,2*n)/(2*m*(2*m-1))
    right=S.binomial(2*n+m-1,m-1)*moment(m,2*n)/2
    assert left==right
print('PASS:64 exact identities equating the two crude rational error bounds')

mp.mp.dps=100
for n in range(1,4):
    mat=S.Matrix([[moment(1,i+j) for j in range(n)] for i in range(n)])
    rhs=S.Matrix([-moment(1,n+i) for i in range(n)])
    sol=mat.inv()*rhs
    p=t**n+sum(sol[i]*t**i for i in range(n))
    Q=S.expand((-1)**n*p.subs(t,-z))
    quot=S.Poly(S.cancel((Q-Q.subs(z,-t))/(z+t)),t)
    P=S.expand(sum(c*moment(1,k[0]) for k,c in quot.terms()))
    R=S.cancel(P/Q)
    for m in range(1,5):
        derivative=S.cancel(S.diff(R.subs(z,a*a),a,2*m-2)/S.factorial(2*m))
        num,den=S.fraction(derivative)
        assert S.Poly(den,a).monic()==S.Poly(Q.subs(z,a*a)**(2*m-1),a).monic()
        assert all(c==0 for k,c in enumerate(reversed(S.Poly(num,a).all_coeffs())) if k%2)
        aa=S.Rational(2*n+1,2)
        value=S.harmonic(n,2*m+1)+derivative.subs(a,aa)
        bound=S.binomial(4*n+2*m-1,2*m-2)*moment(1,2*n)/(2*m*(2*m-1)*aa**(4*n+2*m))
        err=abs(mp.zeta(2*m+1)-mp.mpf(int(value.p))/int(value.q))
        assert err<=mp.mpf(int(bound.p))/int(bound.q)
print('PASS:12 reduced derivative denominatorsQ^(2m−1) and numerical rational-enclosure corroborations')

m=S.symbols('m',positive=True)
f=S.Function('mu')(t)
aa=m+S.Rational(1,2)
nextdensity=((2*m-1)*f-2*t*S.diff(f,t))/(2*m+1)
assert S.simplify(aa*t**(-aa)*nextdensity+S.diff(t**(1-aa)*f,t))==0
print('PASS: generic Beta-product density identity')
