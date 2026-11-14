#!/usr/bin/env python3
"""Exact saddle algebra plus positive-integral diagnostics for the proved rate."""
from pathlib import Path
from math import comb,lcm
import json
import sympy as S
import mpmath as mp

v,u,z=S.symbols('v u z',positive=True)
# Exact tangent numerator identities and root isolation, independent of decimals.
for m,expected in [(1,2*z*(z**4+10*z*z-7)),
                   (2,4*z*(z**6+7*z**4-21*z*z+5))]:
    imaginary=S.expand(S.im(S.expand((z+S.I)**(2*m+3)*(z-3*S.I))))
    assert S.expand(imaginary-expected)==0
assert S.Poly(u**3+7*u*u-21*u+5,u).count_roots(2,3)==1
assert (u**3+7*u*u-21*u+5).subs(u,2)==-1
assert (u**3+7*u*u-21*u+5).subs(u,3)==32
assert S.simplify((S.sqrt(2)-1)**4-(17-12*S.sqrt(2)))==0
assert S.Rational(32)*S.Rational(9,128)==S.Rational(9,4)>1
for n in range(1,101):
    dn=lcm(*range(1,n+1))
    assert all(dn%comb(n-1,k)==0 for k in range(n))
    assert dn*n>=2**(n-1)
print('PASS exact saddle polynomials, cubic root isolation, Apéry value and lcm lower bounds')

mp.mp.dps=55

def saddle(m):
    lo=mp.mpf('1e-40');hi=mp.mpf(m+3)
    for _ in range(230):
        mid=(lo+hi)/2
        f=(2*m+3)*mp.atan(1/(2*mid))-mp.atan(3/(2*mid))-mp.pi
        if f>0:lo=mid
        else:hi=mid
    return (lo+hi)/2

def phi(m,v):
    J=lambda w:w*mp.log(v*v+w*w)-2*w+2*v*mp.atan(w/v)
    return 2*m+(2*m+3)*J(mp.mpf('.5'))-J(mp.mpf('1.5'))-2*mp.pi*v

def density(m,y):
    coeff={1:1}
    for _ in range(1,m):
        nxt={}
        for j,a in coeff.items():
            nxt[j]=nxt.get(j,0)+4*j*j*a
            nxt[j+1]=nxt.get(j+1,0)+2*(j+1)*(2*j+1)*a
        coeff=nxt
    c=mp.csch(mp.pi*y)**2
    return 2*mp.pi**(2*m)/mp.factorial(2*m)*y**(2*m+1)*mp.coth(mp.pi*y)*sum(a*c**j for j,a in coeff.items())

def positive_form(m,n):
    r=n//2
    def integrand(v):
        if v==0:return mp.mpf(0)
        y=n*v;t=y*y
        logarithm=mp.log(t)+(2*m+2)*sum(mp.log(t+k*k) for k in range(1,r))
        logarithm+=(2*m+1)*mp.log(t+r*r)
        logarithm-=sum(mp.log(t+k*k) for k in range(r+1,3*r+2))
        logarithm-=2*m*mp.loggamma(n+1)
        return n*mp.exp(logarithm)*density(m,y)
    vm=saddle(m)
    cuts=[0,vm/2,vm,2*vm,max(4,4*vm),mp.inf]
    return mp.quad(integrand,cuts)

out=[]
for m in range(1,5):
    vm=saddle(m);ph=phi(m,vm)
    rho=(vm*vm+mp.mpf('.25'))**(mp.mpf(2*m+3)/2)/(vm*vm+mp.mpf('2.25'))**mp.mpf('1.5')
    assert abs(mp.log(rho)-ph)<mp.mpf('1e-48')
    if m==1:assert abs(rho-(17-12*mp.sqrt(2)))<mp.mpf('1e-48')
    if m==2:
        uu=4*vm*vm
        assert 2<uu<3 and abs(uu**3+7*uu**2-21*uu+5)<mp.mpf('1e-47')
        assert rho>mp.mpf(9)/128
    diagnostics=[]
    for n in [4,12,30]:
        value=positive_form(m,n)
        assert value>0
        lograte=mp.log(value)/n
        diagnostics.append(dict(n=n,log_positive_form=str(mp.log(value)),finite_log_rate=str(lograte)))
        print('positive integral m,n=',m,n,'log(L)/n=',mp.nstr(lograte,15),'limit=',mp.nstr(ph,15))
    out.append(dict(m=m,saddle=str(vm),proved_log_rate=str(ph),proved_root_rate=str(rho),diagnostics=diagnostics))
Path(__file__).with_name('balanced_diagonal_positive_rate_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS exact saddle algebra and twelve positive-integral diagnostics; convergence is proved in the note, not inferred here')
