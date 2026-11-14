#!/usr/bin/env python3
"""Check generalized Wilson Rodrigues against independent Bernoulli moments.

The mathematical all-index proof is in balanced-generalized-wilson-rodrigues.md.
This checker also tests nonzero contour-adjoint examples numerically, rather
than only comparing vanishing orthogonality conditions.
"""
from math import factorial, comb
from pathlib import Path
import json
import sympy as S
import mpmath as mp

x, t = S.symbols('x t')

def raise_polynomial(f, m, s):
    a = S.Rational(s, 2)
    A = (x+a)*(x+a+1)**(2*m+1)
    numerator = S.Poly(S.expand(A*f.subs(x,x+S.Rational(1,2))
                       -A.subs(x,-x)*f.subs(x,x-S.Rational(1,2))),x)
    assert numerator.nth(0)==0
    q, r = S.div(numerator, S.Poly(2*x,x))
    assert r.is_zero
    assert all(q.nth(k)==0 for k in range(1,q.degree()+1,2))
    return q.as_expr()

def rodrigues(m,n):
    f=S.Integer(1)
    for s in reversed(range(n)):
        f=raise_polynomial(f,m,s)
    D=m*n+m+1
    assert S.Poly(f,x).LC()==S.rf(D,n)
    p=sum(S.Poly(f,x).nth(2*k)*(-t)**k for k in range(m*n+1))
    return f,S.expand((-1)**(m*n)*p/S.rf(D,n))

def moment(j,e):
    return (-1)**e*S.bernoulli(2*e+2,0)*S.binomial(2*e+2*j+1,2*j)/(2*e+2)

rows=[]
for m,n in [(m,n) for m in range(1,5) for n in range(1,5)]+[(5,2),(2,6),(3,5)]:
    f,p=rodrigues(m,n)
    pp=S.Poly(p,t)
    assert pp.degree()==m*n and pp.LC()==1
    residuals=[S.factor(sum(pp.nth(e)*moment(j,m-j+k+e) for e in range(m*n+1)))
               for j in range(1,m+1) for k in range(n)]
    assert all(v==0 for v in residuals)
    assert f.subs(x,0)==factorial(n)**(2*m+1)
    assert f.subs(x,1)==factorial(n)*factorial(n+1)**(2*m)
    expected=S.Rational((-1)**(m*n)*factorial(n)**(2*m),comb((m+1)*n+m,n))
    assert pp.nth(0)==expected
    assert p.subs(t,-1)==expected*(n+1)**(2*m)
    rows.append(dict(m=m,n=n,degree=m*n,orthogonality_checks=m*n,
                     p_coefficients_low_to_high=[str(pp.nth(e)) for e in range(m*n+1)]))
    print('PASS Rodrigues / all rational moments / both evaluations:',m,n)

# Leading-coefficient and endpoint identities on arbitrary polynomials.
for m in range(1,5):
    for s in range(5):
        f=3*x**8-2*x**4+7*x*x-5
        b=raise_polynomial(f,m,s)
        assert S.Poly(b,x).LC()==3*(4+(m+1)*s+2*m+1)
        if s==0:
            assert b.subs(x,0)==f.subs(x,S.Rational(1,2))
        else:
            assert b.subs(x,S.Rational(s,2))==(s+1)**(2*m+1)*f.subs(x,S.Rational(s+1,2))
        assert b.subs(x,S.Rational(s+2,2))==(s+1)*(s+2)**(2*m)*f.subs(x,S.Rational(s+3,2))
print('PASS 20 arbitrary-polynomial leading / endpoint cases')

mp.mp.dps=45

def weight(m,s,z):
    if z==0:
        return mp.mpf(4) if s==0 else mp.mpf(0)
    a=mp.mpf(s)/2
    return mp.gamma(a+z)*mp.gamma(a-z)*(mp.gamma(1+a+z)*mp.gamma(1+a-z))**(2*m+1)/(mp.gamma(2*z)*mp.gamma(-2*z))

for m in range(1,4):
    for s in range(4):
        z=mp.mpc('0.17','0.83')
        A=lambda u:(u+mp.mpf(s)/2)*(u+mp.mpf(s)/2+1)**(2*m+1)
        first=weight(m,s+1,z+mp.mpf('.5'))/weight(m,s,z)
        second=weight(m,s+1,z-mp.mpf('.5'))/weight(m,s,z)
        assert abs(first+(2*z+1)*A(z)/(2*z))<mp.mpf('1e-37')
        assert abs(second+(2*z-1)*A(-z)/(2*z))<mp.mpf('1e-37')
print('PASS 24 gamma Pearson-ratio checks')

# Nonzero adjoint identities: integrate over y>0, dividing the common 2i.
for m,s,j,k in [(1,0,1,2),(2,0,1,2),(2,0,2,1),(2,2,1,2),(3,1,2,1)]:
    f_expr=1+2*x*x
    b_expr=raise_polynomial(f_expr,m,s)
    f=S.lambdify(x,f_expr,'mpmath')
    b=S.lambdify(x,b_expr,'mpmath')
    r=m-j
    g=lambda z:mp.sin(mp.pi*z)**(2*r)*z**(2*k)
    dg=lambda z:(g(z+mp.mpf('.5'))-g(z-mp.mpf('.5')))/(2*z)
    left=mp.quad(lambda y:mp.re(g(1j*y)*weight(m,s,1j*y)*b(1j*y)),[0,1,3,8,mp.inf])
    right=mp.quad(lambda y:mp.re(dg(1j*y)*weight(m,s+1,1j*y)*f(1j*y)),[0,1,3,8,mp.inf])
    rel=abs(left-right)/max(1,abs(left),abs(right))
    assert rel<mp.mpf('1e-32'), (m,s,j,k,left,right,rel)
    assert abs(left)>mp.mpf('1e-8')
    print('PASS nonzero contour-adjoint quadrature:',m,s,j,k,mp.nstr(rel,4))

Path(__file__).with_name('balanced_generalized_wilson_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
print('PASS all finite Rodrigues, exact-moment, endpoint, gamma, and adjoint checks')
