#!/usr/bin/env python3
"""Exact Jacobi normalizations plus independent Euler-kernel quadrature."""
from fractions import Fraction as F
from math import comb
import sympy as s
import mpmath as mp

v=s.symbols('v')
JJ=lambda f:s.cancel(-s.diff(v*(1-v)**2*s.diff(f,v),v)/(1-v))
coeffs=[]
for n in range(31):
    cc=[(-1)**(n+j)*comb(n,j)*comb(n+j+1,j) for j in range(n+1)]
    coeffs.append(cc)
    assert cc[0]==(-1)**n
    norm=sum((F(cc[i]*cc[j],(i+j+1)*(i+j+2))
              for i in range(n+1) for j in range(n+1)),F(0))
    assert norm==F(1,2*(n+1))
    logmoment=sum((F(c,(j+1)**2) for j,c in enumerate(cc)),F(0))
    assert logmoment==F((-1)**n,(n+1)**2)
    if n<=12:
        poly=sum(c*v**j for j,c in enumerate(cc))
        assert s.expand(JJ(poly)-n*(n+2)*poly)==0
        assert s.expand(poly-s.jacobi(n,1,0,2*v-1))==0
for n in range(13):
    for m in range(n):
        val=sum((F(coeffs[n][i]*coeffs[m][j],(i+j+1)*(i+j+2))
                 for i in range(n+1) for j in range(m+1)),F(0))
        assert val==0
print('31 exact endpoints, norms and logarithmic moments; 13 eigen-equations; 78 orthogonal pairs: passed.',flush=True)

# The finite correction to K0 is exactly minus its first n Jacobi modes.
for n in range(1,31):
    direct=[F(0)]*n
    spectral=[F(0)]*n
    for j in range(1,n+1):
        cj=(-1)**j*comb(n,j)*comb(n+j,j)
        for r in range(j):
            direct[r]+=F(cj,j)
    for ell in range(1,n+1):
        for r,cc in enumerate(coeffs[ell-1]):
            spectral[r]-=F(2*(-1)**(ell-1)*cc,ell)
    assert direct==spectral
print('30 exact Legendre-to-Jacobi projection identities: passed.',flush=True)

mp.mp.dps=60
def kernel(q,u):
    return mp.gamma(1-q)*mp.gamma(1+q)*mp.hyp2f1(1-q,1+q,2,1-u)
def poly(n,u):
    return mp.polyval(list(reversed(coeffs[n])),u)
for q in [mp.mpf(1)/3,mp.mpf(1)/2,mp.mpf(3)/4]:
    errors=[]
    for n in range(4):
        integrand=lambda w:mp.exp(-w)*(1-mp.exp(-w))*kernel(q,mp.exp(-w))*poly(n,mp.exp(-w))
        val=mp.quad(integrand,[0,1,5,20,50,100])
        target=(-1)**n/((n+1)**2-q*q)
        err=abs(val-target)
        assert err<mp.mpf('1e-34'),(q,n,err)
        errors.append(err)
    print('Four Jacobi kernel coefficients at tau=',mp.nstr(q,8),
          ' max quadrature error=',mp.nstr(max(errors),6),flush=True)
print('All checks passed.',flush=True)
