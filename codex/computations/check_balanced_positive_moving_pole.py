#!/usr/bin/env python3
"""Check positive moving-pole transport and integer-valued numerator arithmetic."""
from fractions import Fraction as Q
from math import comb, factorial, lcm
from pathlib import Path
import json
import sympy as S
import mpmath as mp

base=Path(__file__).parent
rows=json.loads((base/'balanced_rodrigues_integrality_checks.json').read_text())
cache={}
def bern(k):
    if k not in cache:cache[k]=Q(S.bernoulli(k,0))
    return cache[k]
def moment(m,e):return (-1)**e*bern(2*e+2)*Q(comb(2*e+2*m+1,2*m),2*e+2)
def ev(c,x):return sum(v*x**e for e,v in enumerate(c))
def quotient(c,a):
    out=[Q(0)]*(len(c)-1);out[-1]=c[-1]
    for e in range(len(out)-2,-1,-1):out[e]=c[e+1]+a*out[e+1]
    return out

def associated(row,a):
    m,n=row['m'],row['n'];d=m*n;H=lcm(*range(1,n+1))**(2*m-2)
    q=[Q(c) for c in row['integer_polynomial_coefficients_low_to_high']]
    norm=[v/q[0] for v in q];R=ev(norm,-a*a)
    g=quotient(norm,-a*a)
    assert norm[0]-a*a*g[0]==R
    mu=sum(v*moment(m,e) for e,v in enumerate(g))
    harmonic=sum(Q(1,k**(2*m+1)) for k in range(1,a+1))
    approx=harmonic+Q(1,2*m*a**(2*m))-Q(1,2*a**(2*m+1))-mu/(a**(2*m)*R)
    f=[Q(0)]*(2*d+1)
    for e,v in enumerate(norm):f[2*e]=H*(-1)**e*v
    fa=ev(f,a)
    assert fa.denominator==1
    P=[Q(0)]*(2*m)+quotient(f,a)
    L=2*d+2*m-1;r=2*m-1
    derivative=[P[e+r]*Q(factorial(e+r),factorial(e)) for e in range(L-r+1)]
    Bvalue=sum(c*bern(e) for e,c in enumerate(derivative))
    b_formula=a**(2*m)*fa*harmonic+fa/(2*m)-Q(H,2*a)+Bvalue/factorial(2*m)
    assert b_formula==H*a**(2*m)*R*approx
    d2d=lcm(*range(1,2*d+1));dL=lcm(*range(1,L+1));d2d1=lcm(*range(1,2*d+2))
    assert all((d2d*ev(P,k)).denominator==1 for k in range(L+1))
    assert (d2d*dL**r*d2d1*Bvalue).denominator==1
    J=max(a,L);M=factorial(2*m)*lcm(*range(1,J+1))**(2*m+1)
    assert (M*b_formula).denominator==1
    return approx,norm

checks=0
for row in rows:
    if row['m']>4 or row['n']>5:continue
    for a in sorted(set([1,row['n']//2+1,row['n']+1,2*row['n']+1])):
        associated(row,a);checks+=1
print('PASS',checks,'exact single-quotient/reflection/associated-clearing cases')

mp.mp.dps=60

def density(m,y):
    coeff={1:1}
    for u in range(1,m):
        nxt={}
        for j,v in coeff.items():
            nxt[j]=nxt.get(j,0)+4*j*j*v
            nxt[j+1]=nxt.get(j+1,0)+2*(j+1)*(2*j+1)*v
        coeff=nxt
    c=1/mp.sinh(mp.pi*y)**2
    return 2*mp.pi**(2*m)/factorial(2*m)*y**(2*m+1)*mp.coth(mp.pi*y)*sum(v*c**j for j,v in coeff.items())

quads=0
for m,n,a in [(1,2,2),(2,2,2),(2,2,4),(2,4,3),(3,2,2)]:
    row=next(v for v in rows if v['m']==m and v['n']==n)
    approx,norm=associated(row,a);r=n//2
    q0=Q(row['integer_polynomial_coefficients_low_to_high'][0],row['binomial_normalizer'])
    p=[q0*v for v in norm]
    pmp=[mp.mpf(v.numerator)/v.denominator for v in p]
    pe=lambda t:sum(v*t**e for e,v in enumerate(pmp))
    Sp=lambda t:t*mp.fprod((t+k*k)**(2*m+2) for k in range(1,r))*(t+r*r)**(2*m+1)
    Db=lambda t:mp.fprod(t+(a+k)**2 for k in range(-r,r+1))
    left=mp.quad(lambda y:pe(y*y)/(y*y+a*a)*density(m,y),[0,1,3,8,mp.inf])
    right=mp.factorial(n)/mp.rf(m*n+m+1,n)*mp.quad(lambda y:Sp(y*y)/Db(y*y)*density(m,y),[0,1,3,8,mp.inf])
    assert left>0 and right>0 and abs(left-right)<mp.mpf('1e-45')*max(1,abs(left),abs(right))
    error=mp.zeta(2*m+1)-mp.mpf(approx.numerator)/approx.denominator
    predicted=left/(a**(2*m)*pe(-a*a))
    assert error>0 and abs(error-predicted)<mp.mpf('1e-45')*max(1,abs(error))
    print('PASS positive transport and zeta remainder:',m,n,a,mp.nstr(error,12))
    quads+=1

(base/'balanced_positive_moving_pole_checks.json').write_text(json.dumps(dict(
    exact_associated_arithmetic_cases=checks,positive_quadratures=quads),indent=2)+'\n')
print('PASS all moving-pole transport and associated arithmetic checks')
