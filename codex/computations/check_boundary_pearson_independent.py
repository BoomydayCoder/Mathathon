#!/usr/bin/env python3
"""Independent direct real-integral check of the Wilson boundary normalization.
Uses SymPy to differentiate polynomials and mpmath for quadrature, not the
rational moment or partial-fraction routines in the author's checker.
"""
import sympy as S
import mpmath as mp
mp.mp.dps=75
x,y=S.symbols('x y'); count=0

def op(p,M):
    a=(x+1)**3*(x-M)/(2*x*(2*x+1))
    return S.cancel(a*(p.subs(x,x+1)-p)+a.subs(x,-x)*(p.subs(x,x-1)-p))
def integrand(p,M):
    real=S.Poly(S.expand(p),x)
    expr=sum(c*(-y*y)**(k//2) for (k,),c in real.terms())
    f=S.lambdify(y,expr,'mpmath')
    def fun(v):
        density=1/mp.pi if not v else mp.pi**2*v**3*mp.coth(mp.pi*v)/mp.sinh(mp.pi*v)**2
        return f(v)*density/mp.fprod(v*v+j*j for j in range(1,M+1))
    return fun
for M in [1,2,3,4,6,10]:
    p=x*x+S.Rational(2,3);q=x**6-S.Rational(5,7)
    v=mp.quad(integrand(S.expand(p*op(q,M)-q*op(p,M)),M),[0,mp.mpf('0.5'),2,8,mp.inf])
    c=S.Rational(M,2*S.factorial(M)**2)
    expected=-c*(p.subs(x,0)*q.subs(x,1)-p.subs(x,1)*q.subs(x,0))
    e=mp.mpf(str(S.N(expected,80)))
    assert abs(v-e)<mp.mpf('1e-60')*max(abs(e),mp.mpf('1e-30'))
    count+=1
    for k in range(max(0,M-2)):
        Q=S.prod(j*j-x*x for j in range(1,k+1))/(S.factorial(k)*S.factorial(k+1))
        v=mp.quad(integrand(Q,M),[0,mp.mpf('0.5'),2,8,mp.inf])
        expected=c*(M-k-2-M*(k+1)*(S.harmonic(M-1,2)-S.harmonic(k+1,2)))/((M-k-2)*(M-k-1))
        e=mp.mpf(str(S.N(expected,80)))
        assert abs(v-e)<mp.mpf('1e-60')*max(abs(e),mp.mpf('1e-30'))
        count+=1
print(f'PASS: {count} direct positive-axis quadratures reproduce the boundary defect and closed Newton moments to relative error below 1e-60')
